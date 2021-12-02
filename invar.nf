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
}
