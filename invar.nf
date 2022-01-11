#!/usr/bin/env nextflow

/*
 * Main alignment work flow.
 */

nextflow.enable.dsl = 2

include { invar12 } from './processes/1_parse/invar12'
include { invar34 } from './processes/1_parse/invar34'
include { sizeAnnotation } from './processes/2_size_annotation/sizeAnnotation'

/*
 * Mini work flow for part one (parsing).
 */
workflow parse
{
    take:
        bamChannel
        tumourMutationsChannel
        layoutChannel

    main:
        invar12(bamChannel, tumourMutationsChannel)
        invar34(invar12.out, tumourMutationsChannel, layoutChannel)

    emit:
        mutationsFile = invar34.out
}

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

    sizeAnnotation(bamChannel, parse.out, tumourMutationsChannel)
}
