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
include { runAnalysis } from '../processes/5_analysis'

include { rdsDiff as diff1; rdsDiff as diff2; rdsDiff as diff3; rdsDiff as diff4; rdsDiff as diff5; rdsDiff as diff6;
          rdsDiff as diff7; rdsDiff as diff8; rdsDiff as diff9; rdsDiff as diff10; diff as diff11 } from '../testing/diff'

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
        tuple pname, it, file("testdata/${pname}/reference/REFERENCE_${it.baseName}.?sv", checkIfExists: true)
    }
}

workflow
{
    // dumpParams(log, params)

    tumourMutationsChannel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
    layoutChannel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

    // createMutationsTable

    createMutationsTable(channel.fromPath('testdata/createMutationsTable/source/mutation_table.tsv'),
                         tumourMutationsChannel)

    mapForDiff('createMutationsTable', createMutationsTable.out.filteredMutationsFile) | diff1

    // offTargetErrorRates

    offTargetErrorRates(channel.fromPath("testdata/offTargetErrorRates/source/mutation_table.filtered.rds"),
                        layoutChannel)

    offTargetOutputs =
        offTargetErrorRates.out.locusErrorRates
            .mix(offTargetErrorRates.out.cosmicErrorRates)
            .mix(offTargetErrorRates.out.noCosmicErrorRates)

    mapForDiff('offTargetErrorRates', offTargetOutputs) | diff2

    // createOnTargetMutationsTable

    createOnTargetMutationsTable(channel.fromPath('testdata/createOnTargetMutationsTable/source/mutation_table.filtered.rds'),
                                 tumourMutationsChannel,
                                 layoutChannel,
                                 channel.fromPath('testdata/createOnTargetMutationsTable/source/error_rates.off_target.no_cosmic.rds'))

    mapForDiff('createOnTargetMutationsTable',
        createOnTargetMutationsTable.out.onTargetMutationsFile.mix(createOnTargetMutationsTable.out.backgroundErrorRates)) | diff3

    // onTargetErrorRatesAndFilter

    onTargetErrorRatesAndFilter(channel.fromPath('testdata/onTargetErrorRatesAndFilter/source/mutation_table.on_target.all.rds'),
                                tumourMutationsChannel)

    mapForDiff('onTargetErrorRatesAndFilter', onTargetErrorRatesAndFilter.out.locusErrorRates.mix(onTargetErrorRatesAndFilter.out.onTargetMutationsFile)) | diff4

    // sizeAnnotation

    sizeAnnotationInsertsChannel = channel.of(['SLX-19721:SXTLI001']).combine(
        channel.fromPath("testdata/annotateMutationsWithFragmentSize/source/SLX19721SXTLI001.inserts.tsv"))

    annotateMutationsWithFragmentSize(sizeAnnotationInsertsChannel,
                                      channel.fromPath('testdata/annotateMutationsWithFragmentSize/source/mutation_table.on_target.rds'))

    annotateMutationsWithFragmentSize_perPatientChannel =
        annotateMutationsWithFragmentSize.out.mutationsFiles
            .flatMap {
                sampleId, indexFile, mutationsFiles ->
                combineIndexFileWithMutationsFiles(indexFile, mutationsFiles)
            }
            .filter {
                sampleId, pt, f ->
                sampleId == 'SLX-19721:SXTLI001' && pt in ['PARA_002', 'PARA_028']
            }
            .map { sampleId, pt, f -> f }

    mapForDiff('annotateMutationsWithFragmentSize', annotateMutationsWithFragmentSize_perPatientChannel) | diff5

    // Outlier suppression

    markOutliersChannel = channel.of(['SLX-19721:SXTLI001', 'PARA_002']).combine(
        channel.fromPath("testdata/markOutliers/source/mutation_table.with_sizes.SLX19721SLXLI001.PARA_002.rds"))

    markOutliers(markOutliersChannel)

    mapForDiff('markOutliers', markOutliers.out.mutationsFile.map { s, pt, f -> f }) | diff6

    // Size characterisation

    sizeCharacterisation(channel.fromPath("testdata/sizeCharacterisation/source/*.rds").collect())

    mapForDiff('sizeCharacterisation', sizeCharacterisation.out.allSizesFile.mix(sizeCharacterisation.out.summaryFile)) | diff7

    // Annotate main mutations with outlier suppression flags.

    annotateMutationsWithOutlierSuppression(
        channel.fromPath("testdata/annotateMutationsWithOutlierSuppression/source/mutation_table.on_target.rds"),
        channel.fromPath("testdata/annotateMutationsWithOutlierSuppression/source/mutation_table.outliersuppressed.*.rds").collect())

    mapForDiff('annotateMutationsWithOutlierSuppression', annotateMutationsWithOutlierSuppression.out.mutationsFile) | diff8

    // Generalised Likelihood Ratio Test

    generalisedLikelihoodRatioTestChannelSpecific = channel.of(['SLX-19721:SXTLI001', 'PARA_002']).combine(
        channel.fromPath("testdata/generalisedLikelihoodRatioTest/source/mutation_table.outliersuppressed.SLX19721SXTLI001.1.rds"))

    generalisedLikelihoodRatioTestChannelNonSpecific = channel.of(['SLX-19721:SXTLI001', 'PARA_028']).combine(
        channel.fromPath("testdata/generalisedLikelihoodRatioTest/source/mutation_table.outliersuppressed.SLX19721SXTLI001.3.rds"))

    generalisedLikelihoodRatioTestSizeChannel =
        channel.fromPath("testdata/generalisedLikelihoodRatioTest/source/size_characterisation.rds")

    generalisedLikelihoodRatioTestSpecific(
        generalisedLikelihoodRatioTestChannelSpecific, generalisedLikelihoodRatioTestSizeChannel)

    generalisedLikelihoodRatioTestNonSpecific(
        generalisedLikelihoodRatioTestChannelNonSpecific, generalisedLikelihoodRatioTestSizeChannel)

    mapForDiff('generalisedLikelihoodRatioTest', generalisedLikelihoodRatioTestSpecific.out.invarScores.map { s, pt, f -> f }) | diff9
    mapForDiff('generalisedLikelihoodRatioTest', generalisedLikelihoodRatioTestNonSpecific.out.invarScores.map { s, pt, f -> f }) | diff10

    // Mutation tracking (part of analysis)

    runAnalysis(channel.fromPath("testdata/runAnalysis/source/mutation_table.rds"),
                tumourMutationsChannel,
                layoutChannel,
                channel.fromPath("testdata/runAnalysis/source/background_error_rates.rds"),
                channel.fromPath("testdata/runAnalysis/source/locus_error_rates.on_target.rds"),
                channel.fromPath("testdata/runAnalysis/source/error_rates.off_target.no_cosmic.rds"),
                channel.fromPath("testdata/runAnalysis/source/size_characterisation.rds"),
                channel.fromPath("testdata/runAnalysis/source/invar_scores.rds"))

    mapForDiff('runAnalysis', runAnalysis.out.mutationsTracking) | diff11
}
