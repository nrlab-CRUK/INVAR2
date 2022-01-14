#!/usr/bin/env nextflow

/*
 * INVAR testing work flow.
 */

nextflow.enable.dsl = 2

include { createMutationsTable } from './processes/1_parse'

def compareFiles(logger, process, generated, reference)
{
    generated.withReader
    {
        greader ->
        
        reference.withReader
        {
            rreader ->
            
            logger.warn "Compare ${process} output to reference."
    
            def line = 0
            
            while (true)
            {
                gline = greader.readLine()
                rline = rreader.readLine()
                ++line
            
                if (!gline && !rline)
                {
                    logger.warn "${process} files are the same"
                    break
                }
                if (!gline && rline)
                {
                    logger.error "${process}: there are fewer lines in the reference than the generated file."
                    logger.error "${process} reference ${reference}"
                    logger.error "${process} generated ${generated.value}"
                    break
                }
                if (gline && !rline)
                {
                    logger.error "${process}: there are more lines in the reference than the generated file."
                    logger.error "${process} reference ${reference}"
                    logger.error "${process} generated ${generated.value}"
                    break
                }
                if (gline != rline)
                {
                    logger.error "${process}: files differ on line ${line}:"
                    logger.error rline
                    logger.error gline
                    logger.error "${process} reference ${reference}"
                    logger.error "${process} generated ${generated.value}"
                    break
                }
            }
        }
    }
}

workflow
{
    tumourMutationsChannel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
    layoutChannel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)
    
    createMutationsTable(channel.fromPath('testing/createMutationsTable/source/mutation_table.tsv'),
                         tumourMutationsChannel,
                         layoutChannel)

    compareFiles(log, 'createMutationsTable',
                 createMutationsTable.out.filteredMutationsTSV.first(),
                 file('testing/createMutationsTable/reference/createMutationsTable.out.tsv', checkIfExists: true))
    
}
