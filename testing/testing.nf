#!/usr/bin/env nextflow

/*
 * INVAR testing work flow.
 */

nextflow.enable.dsl = 2

include { createMutationsTable; offTargetErrorRates;
          createOnTargetMutationsTable; onTargetErrorRatesAndFilter } from '../processes/1_parse'
include { annotateMutationsWithFragmentSize  } from '../processes/2_size_annotation'
include { markOutliers  } from '../processes/3_outlier_suppression'

include { diff as diff1; diff as diff2; diff as diff3; diff as diff4; diff as diff5; diff as diff6 } from './diff'

def dumpParams(logger, params)
{
    def keys = params.keySet().sort()
    for (k in keys)
    {
        logger.warn "${k} = '${params[k]}'"
    }
}

def compareFiles(logger, process, generated, reference)
{
    generated.withReader
    {
        greader ->

        reference.withReader
        {
            rreader ->

            def line = 0

            while (true)
            {
                gline = greader.readLine()
                rline = rreader.readLine()
                ++line

                if (!gline && !rline)
                {
                    logger.info "${process} ${generated.name}: files are the same."
                    break
                }
                if (!gline && rline)
                {
                    logger.error "${process} ${generated.name}: there are fewer lines in the reference than the generated file.\n${reference}\n${generated}"
                    break
                }
                if (gline && !rline)
                {
                    logger.error "${process} ${generated.name}: there are more lines in the reference than the generated file.\n${reference}\n${generated}"
                    break
                }
                if (gline != rline)
                {
                    logger.error "${process} ${generated.name}: files differ on line ${line}:\n${rline}\n${gline}\n${reference}\n${generated}"
                    break
                }
            }
        }
    }
}

def mapForDiff(pname, channel)
{
    channel.map
    {
        tuple pname, it, file("testdata/${pname}/reference/REFERENCE_${it.name}", checkIfExists: true)
    }
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

    mapForDiff('offTargetErrorRates', offTargetErrorRates.out.locusErrorRatesTSV.mix(offTargetErrorRates.out.errorRatesTSV)) | diff2

    // createOnTargetMutationsTable

    createOnTargetMutationsTable(channel.fromPath('testdata/createOnTargetMutationsTable/source/mutation_table.filtered.rds'),
                                 tumourMutationsChannel,
                                 layoutChannel,
                                 channel.fromPath('testdata/createOnTargetMutationsTable/source/mutation_table.error_rates.no_cosmic.rds'))

    mapForDiff('createOnTargetMutationsTable', createOnTargetMutationsTable.out.onTargetMutationsTSV) | diff3

    // onTargetErrorRatesAndFilter

    onTargetErrorRatesAndFilter(channel.fromPath('testdata/onTargetErrorRatesAndFilter/source/mutation_table.on_target.all.rds'),
                                tumourMutationsChannel,
                                layoutChannel)

    mapForDiff('onTargetErrorRatesAndFilter', onTargetErrorRatesAndFilter.out.locusErrorRatesTSV.mix(onTargetErrorRatesAndFilter.out.onTargetMutationsTSV)) | diff4

    // sizeAnnotation

    sizeAnnotationInsertsChannel = channel.of(['SLX-19721', 'SXTLI001']).combine(
        channel.fromPath("testdata/annotateMutationsWithFragmentSize/source/SLX-19721_SXTLI001.inserts.tsv"))

    annotateMutationsWithFragmentSize(sizeAnnotationInsertsChannel,
                                      channel.fromPath('testdata/annotateMutationsWithFragmentSize/source/mutation_table.on_target.rds'))

    mapForDiff('annotateMutationsWithFragmentSize', annotateMutationsWithFragmentSize.out.mutationsTSV) | diff5

    // Outlier suppression

    markOutliersChannel = channel.of(['SLX-19721', 'SXTLI001']).combine(
        channel.fromPath("testdata/markOutliers/source/combined.polished.size_ann.SLX-19721_SXTLI001.rds"))

    markOutliers(markOutliersChannel)

    mapForDiff('markOutliers', markOutliers.out.mutationsTSV) | diff6
}
