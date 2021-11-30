#!/usr/bin/env nextflow

/*
 * Main alignment work flow.
 */

nextflow.enable.dsl = 2

include { invar1 } from './processes/invar1'

/*
 * Main work flow.
 */
workflow
{
    invar1()
    
    // This next bit is all INVAR2.sh does.
    
    all_mutations =
        invar1.out.mutationFiles.collectFile(name: "${params.FINAL_PREFIX}.combined.final.ann.tsv",
                                             storeDir: 'mutations', keepHeader: true, skip: 1, sort: false)
}
