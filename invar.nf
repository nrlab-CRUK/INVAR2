#!/usr/bin/env nextflow

/*
 * Main alignment work flow.
 */

nextflow.enable.dsl = 2

include { parse } from './processes/1_parse'
include { sizeAnnotation } from './processes/2_size_annotation'
include { outlierSuppression } from './processes/3_outlier_suppression'

/*
 * Main work flow.
 */
workflow
{
    tumourMutationsChannel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
    layoutChannel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

    bamChannel = channel.fromPath(params.INPUT_FILES, checkIfExists: true)
        .splitCsv(header: true, by: 1, strip: true)
        .map {
            row ->
            bam = file(row.FILE_NAME, checkIfExists: true)
            index = file("${bam.name}.bai") // The index file may or may not exist.
            tuple row.POOL, row.BARCODE, bam, index
        }

    parse(bamChannel, tumourMutationsChannel, layoutChannel)

    sizeAnnotation(bamChannel, parse.out.onTargetMutationsFile, tumourMutationsChannel)

    outlierSuppression(parse.out.onTargetMutationsFile, sizeAnnotation.out.mutationsFiles)
}
