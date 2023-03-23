#!/usr/bin/env nextflow

@Grab('org.apache.commons:commons-lang3:3.12.0')

import static org.apache.commons.lang3.StringUtils.isNotBlank
import static org.apache.commons.lang3.StringUtils.trimToEmpty

/*
 * Main alignment work flow.
 */

nextflow.enable.dsl = 2

include { findFileOnPath } from './functions/paths'
include { validateParameters; validateReferenceFiles } from './functions/validation'
include { parse } from './processes/1_parse'
include { sizeAnnotation } from './processes/2_size_annotation'
include { outlierSuppression } from './processes/3_outlier_suppression'
include { detection } from './processes/4_detection'
include { analysis } from './processes/5_analysis'

/*
 * Check the pipeline is set up without basic errors.
 */
if (!validateParameters(params))
{
    exit 1
}

/*
 * Check the layout file has the required columns and has files for each
 * sample in the study. Also check the tumour mutations file has the
 * required columns.
 */
if (!validateReferenceFiles(params))
{
    exit 1
}

/*
 * Main work flow.
 */
workflow
{
    tumourMutationsChannel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
    layoutChannel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

    bamChannel = layoutChannel
        .splitCsv(header: true, by: 1, strip: true, quote: '"', charset = params.LAYOUT_TABLE_ENCODING)
        .filter {
            row ->
            row.STUDY == params.STUDY &&
            isNotBlank(row.BAM_FILE) &&
            trimToEmpty(row.ACTIVE).toLowerCase() in [ '', 'yes', 'y', 'true', 't' ]
        }
        .map {
            row ->
            bam = file(findFileOnPath(params.BAM_PATH, row.BAM_FILE), checkIfExists: true)
            index = file("${bam}.bai") // The index file may or may not exist.
            tuple row.SAMPLE_ID, bam, index
        }

    parse(bamChannel, tumourMutationsChannel, layoutChannel)

    sizeAnnotation(bamChannel, parse.out.onTargetMutationsFile, tumourMutationsChannel)

    outlierSuppression(parse.out.onTargetMutationsFile, sizeAnnotation.out.mutationsFiles)

    detection(outlierSuppression.out.perSampleMutationsFiles, outlierSuppression.out.sizeCharacterisationFile)

    analysis(outlierSuppression.out.mutationsFile,
             tumourMutationsChannel,
             layoutChannel,
             parse.out.backgroundErrorRatesFile,
             parse.out.onTargetLocusErrorRatesFile,
             parse.out.offTargetErrorRatesNoCosmic,
             outlierSuppression.out.sizeCharacterisationFile,
             detection.out.invarScores)
}
