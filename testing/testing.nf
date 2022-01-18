#!/usr/bin/env nextflow

/*
 * INVAR testing work flow.
 */

nextflow.enable.dsl = 2

include { createMutationsTable; offTargetErrorRates;
          createOnTargetMutationsTable; onTargetErrorRatesAndFilter } from '../processes/1_parse'
include { annotateMutationsWithFragmentSize  } from '../processes/2_size_annotation'

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


workflow
{
    // dumpParams(log, params)

    tumourMutationsChannel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
    layoutChannel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

    // createMutationsTable

    createMutationsTable(channel.fromPath('testdata/createMutationsTable/source/mutation_table.tsv'),
                         tumourMutationsChannel,
                         layoutChannel)

    createMutationsTable.out.filteredMutationsTSV.first()
        .subscribe onNext:
        {
            genFile ->
            refFile = file("testdata/createMutationsTable/reference/${genFile.name}", checkIfExists: true)
            compareFiles(log, "createMutationsTable", genFile, refFile)
        }

    // offTargetErrorRates

    offTargetErrorRates(channel.fromPath("testdata/offTargetErrorRates/source/mutation_table.filtered.rds"),
                        layoutChannel)

    offTargetErrorRates.out.locusErrorRatesTSV.first()
        .subscribe onNext:
        {
            genFile ->
            refFile = file("testdata/offTargetErrorRates/reference/${genFile.name}", checkIfExists: true)
            compareFiles(log, "offTargetErrorRates", genFile, refFile)
        }

    offTargetErrorRates.out.errorRatesTSV.first().flatten()
        .subscribe onNext:
        {
            genFile ->
            refFile = file("testdata/offTargetErrorRates/reference/${genFile.name}", checkIfExists: true)
            compareFiles(log, "offTargetErrorRates", genFile, refFile)
        }

    // createOnTargetMutationsTable

    createOnTargetMutationsTable(channel.fromPath('testdata/createOnTargetMutationsTable/source/mutation_table.filtered.rds'),
                                 tumourMutationsChannel,
                                 layoutChannel,
                                 channel.fromPath('testdata/createOnTargetMutationsTable/source/mutation_table.error_rates.no_cosmic.rds'))

    createOnTargetMutationsTable.out.onTargetMutationsTSV.first()
        .subscribe onNext:
        {
            genFile ->
            refFile = file("testdata/createOnTargetMutationsTable/reference/${genFile.name}", checkIfExists: true)
            compareFiles(log, "createOnTargetMutationsTable", genFile, refFile)
        }

    // onTargetErrorRatesAndFilter

    onTargetErrorRatesAndFilter(channel.fromPath('testdata/onTargetErrorRatesAndFilter/source/mutation_table.on_target.all.rds'),
                                tumourMutationsChannel,
                                layoutChannel)

    onTargetErrorRatesAndFilter.out.locusErrorRatesTSV.first()
        .subscribe onNext:
        {
            genFile ->
            refFile = file("testdata/onTargetErrorRatesAndFilter/reference/${genFile.name}", checkIfExists: true)
            compareFiles(log, "onTargetErrorRatesAndFilter", genFile, refFile)
        }

    onTargetErrorRatesAndFilter.out.onTargetMutationsTSV.first()
        .subscribe onNext:
        {
            genFile ->
            refFile = file("testdata/onTargetErrorRatesAndFilter/reference/${genFile.name}", checkIfExists: true)
            compareFiles(log, "onTargetErrorRatesAndFilter", genFile, refFile)
        }
    
    // sizeAnnotation
    
    sizeAnnotationInsertsChannel = channel.of(['SLX-19721', 'SXTLI001']).combine(
        channel.fromPath("testdata/annotateMutationsWithFragmentSize/source/SLX-19721_SXTLI001.inserts.tsv"))
    
    annotateMutationsWithFragmentSize(sizeAnnotationInsertsChannel,
                                      channel.fromPath('testdata/annotateMutationsWithFragmentSize/source/mutation_table.on_target.rds'))
    
    annotateMutationsWithFragmentSize.out.mutationsTSV.first()
        .subscribe onNext:
        {
            genFile ->
            refFile = file("testdata/annotateMutationsWithFragmentSize/reference/combined.polished.size_ann.SLX-19721_SXTLI001.tsv", checkIfExists: true)
            compareFiles(log, "annotateMutationsWithFragmentSize", genFile, refFile)
        }
}
