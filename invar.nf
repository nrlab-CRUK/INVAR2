#!/usr/bin/env nextflow

/*
 * Main alignment work flow.
 */

nextflow.enable.dsl = 2

include { invar1 } from './processes/invar1'
include { invar3 } from './processes/invar3'

/*
 * Main work flow.
 */
workflow
{
    invar1() | invar3
}
