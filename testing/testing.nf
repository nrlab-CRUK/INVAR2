#!/usr/bin/env nextflow

/*
 * INVAR testing work flow.
 */

nextflow.enable.dsl = 2

include { createMutationsTable; offTargetErrorRates;
          createOnTargetMutationsTable; onTargetErrorRatesAndFilter } from '../processes/1_parse'
include { annotateMutationsWithFragmentSize ; combineIndexFileWithMutationsFiles } from '../processes/2_size_annotation'
include { markOutliers; sizeCharacterisation; annotateMutationsWithOutlierSuppression } from '../processes/3_outlier_suppression'
include { generalisedLikelihoodRatioTest as generalisedLikelihoodRatioTestSpecific;
          generalisedLikelihoodRatioTest as generalisedLikelihoodRatioTestNonSpecific } from '../processes/4_detection'

include { diff as diff1; diff as diff2; diff as diff3; diff as diff4; diff as diff5; diff as diff6;
          diff as diff7; diff as diff8; diff as diff9; diff as diff10 } from './diff'

def dumpParams(logger, params)
{
    def keys = params.keySet().sort()
    for (k in keys)
    {
        logger.warn "${k} = '${params[k]}'"
    }
}

def mapForDiff(pname, channel)
{
    channel.map
    {
        tuple pname, it, file("testdata/${pname}/reference/REFERENCE_${it.name}", checkIfExists: true)
    }
}

process trimGLRT
{
    executor 'local'
    time '2m'
    cpus 1
    memory '64m'

    input:
        path glrtOutputRDS

    output:
        path outputFilename, emit: "trimmedGLRT"

    shell:
        outputFilename = "${glrtOutputRDS.baseName}.trimmed.tsv"

        """
        Rscript --vanilla "!{params.projectHome}/testing/R/glrtTrim.R" \
            "!{glrtOutputRDS}" "!{outputFilename}"
        """
}

workflow
{
    // dumpParams(log, params)

    tumourMutationsChannel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
    layoutChannel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

    // createMutationsTable

    createMutationsTable(channel.fromPath('testdata/createMutationsTable/source/mutation_table.tsv'),
                         tumourMutationsChannel,
                         layoutChannel)

    mapForDiff('createMutationsTable', createMutationsTable.out.filteredMutationsTSV) | diff1

    // offTargetErrorRates

    offTargetErrorRates(channel.fromPath("testdata/offTargetErrorRates/source/mutation_table.filtered.rds"),
                        layoutChannel)

    mapForDiff('offTargetErrorRates', offTargetErrorRates.out.locusErrorRatesTSV.mix(offTargetErrorRates.out.errorRatesTSV.flatten())) | diff2

    // createOnTargetMutationsTable

    createOnTargetMutationsTable(channel.fromPath('testdata/createOnTargetMutationsTable/source/mutation_table.filtered.rds'),
                                 tumourMutationsChannel,
                                 layoutChannel,
                                 channel.fromPath('testdata/createOnTargetMutationsTable/source/mutation_table.error_rates.no_cosmic.rds'))

    mapForDiff('createOnTargetMutationsTable',
        createOnTargetMutationsTable.out.onTargetMutationsTSV.mix(createOnTargetMutationsTable.out.backgroundErrorRatesTSV)) | diff3

    // onTargetErrorRatesAndFilter

    onTargetErrorRatesAndFilter(channel.fromPath('testdata/onTargetErrorRatesAndFilter/source/mutation_table.on_target.all.rds'),
                                tumourMutationsChannel,
                                layoutChannel)

    mapForDiff('onTargetErrorRatesAndFilter', onTargetErrorRatesAndFilter.out.locusErrorRatesTSV.mix(onTargetErrorRatesAndFilter.out.onTargetMutationsTSV)) | diff4

    // sizeAnnotation

    sizeAnnotationInsertsChannel = channel.of(['SLX-19721', 'SXTLI001']).combine(
        channel.fromPath("testdata/annotateMutationsWithFragmentSize/source/SLX-19721.SXTLI001.inserts.tsv"))

    annotateMutationsWithFragmentSize(sizeAnnotationInsertsChannel,
                                      channel.fromPath('testdata/annotateMutationsWithFragmentSize/source/mutation_table.on_target.rds'))

    annotateMutationsWithFragmentSize_perPatientChannel =
        annotateMutationsWithFragmentSize.out.mutationsFiles
            .flatMap {
                pool, barcode, indexFile, mutationsFiles ->
                combineIndexFileWithMutationsFiles(indexFile, mutationsFiles)
            }
            .filter {
                p, b, pt, f ->
                p == 'SLX-19721' && b == 'SXTLI001' && pt in ['PARA_002', 'PARA_028']
            }
            .map { p, b, pt, f -> f }

    // TODO: Produces only RDS files. Need to change to compare TSV.
    // mapForDiff('annotateMutationsWithFragmentSize', annotateMutationsWithFragmentSize_perPatientChannel) | diff5

    // Outlier suppression

    markOutliersChannel = channel.of(['SLX-19721', 'SXTLI001', 'PARA_002']).combine(
        channel.fromPath("testdata/markOutliers/source/mutation_table.with_sizes.SLX-19721.SLXLI001.PARA_002.rds"))

    markOutliers(markOutliersChannel)

    mapForDiff('markOutliers', markOutliers.out.mutationsTSV) | diff6

    // Size characterisation

    sizeCharacterisation(channel.fromPath("testdata/sizeCharacterisation/source/*.rds").collect())

    mapForDiff('sizeCharacterisation', sizeCharacterisation.out.allSizesTSV.mix(sizeCharacterisation.out.summaryTSV)) | diff7

    // Annotate main mutations with outlier suppression flags.

    annotateMutationsWithOutlierSuppression(
        channel.fromPath("testdata/annotateMutationsWithOutlierSuppression/source/mutation_table.on_target.rds"),
        channel.fromPath("testdata/annotateMutationsWithOutlierSuppression/source/mutation_table.outliersuppressed.*.rds").collect())

    mapForDiff('annotateMutationsWithOutlierSuppression', annotateMutationsWithOutlierSuppression.out.mutationsTSV) | diff8

    // Generalised Likelihood Ratio Test

    generalisedLikelihoodRatioTestChannelSpecific = channel.of(['SLX-19721', 'SXTLI001', 'PARA_002']).combine(
        channel.fromPath("testdata/generalisedLikelihoodRatioTest/source/mutation_table.outliersuppressed.SLX-19721.SXTLI001.1.rds"))

    generalisedLikelihoodRatioTestChannelNonSpecific = channel.of(['SLX-19721', 'SXTLI001', 'PARA_028']).combine(
        channel.fromPath("testdata/generalisedLikelihoodRatioTest/source/mutation_table.outliersuppressed.SLX-19721.SXTLI001.3.rds"))

    generalisedLikelihoodRatioTestSizeChannel =
        channel.fromPath("testdata/generalisedLikelihoodRatioTest/source/size_characterisation.rds")

    generalisedLikelihoodRatioTestSpecific(
        generalisedLikelihoodRatioTestChannelSpecific, generalisedLikelihoodRatioTestSizeChannel)

    generalisedLikelihoodRatioTestNonSpecific(
        generalisedLikelihoodRatioTestChannelNonSpecific, generalisedLikelihoodRatioTestSizeChannel)

    // Need to trim the non-specific file to only include the first iteration.
    trimGLRT(generalisedLikelihoodRatioTestNonSpecific.out.invarScoresFile.map { p, b, pt, f -> f } )

    mapForDiff('generalisedLikelihoodRatioTest', generalisedLikelihoodRatioTestSpecific.out.invarScoresTSV) | diff9
    mapForDiff('generalisedLikelihoodRatioTest', trimGLRT.out.trimmedGLRT) | diff10
}
