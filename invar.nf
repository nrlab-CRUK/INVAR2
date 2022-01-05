#!/usr/bin/env nextflow

/*
 * Main alignment work flow.
 */

nextflow.enable.dsl = 2

include { invar12 } from './processes/1_parse/invar12'
include { invar34 } from './processes/1_parse/invar34'

/*
 * Mini work flow for part one (parsing).
 */
workflow parse
{
    main:
        invar12() | invar34

    emit:
        invar34.out
}

/*
 * Main work flow.
 */
workflow
{
    parse()
}
