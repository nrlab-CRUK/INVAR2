process runAnalysis
{
    memory '2g'
    time   '10m'

    publishDir params.ANALYSIS_DIR, mode: 'link', overwrite: true

    input:
        path osMutationsFile
        path layoutFile
        path onTargetErrorRatesFile
        path offTargetErrorRatesFile
        path sizeCharacterisationFile
        path invarScoresFile

    output:
        path "*.pdf", emit: "plots"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/5_analysis/analysis.R" \
            --mutations="!{osMutationsFile}" \
            --layout="!{layoutFile}" \
            --error-rates="!{onTargetErrorRatesFile}" \
            --off-target-error-rates="!{offTargetErrorRatesFile}" \
            --size-characterisation="!{sizeCharacterisationFile}" \
            --invar-scores="!{invarScoresFile}" \
            --study="!{params.STUDY}" \
            --error-suppression="!{params.ERROR_SUPPRESSION_NAME}" \
            --family-size="!{params.FAMILY_SIZE}" \
            --outlier-suppression=!{params.OUTLIER_SUPPRESSION_THRESHOLD}
        """
}


process generateReport
{
    memory '2g'
    time   '10m'

    publishDir params.ANALYSIS_DIR, mode: 'link', overwrite: true

    input:
        path osMutationsFile
        path layoutFile
        path onTargetErrorRatesFile
        path offTargetErrorRatesFile
        path sizeCharacterisationFile
        path invarScoresFile

    output:
        path "*.html"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/5_analysis/knitAnalysis.R" \
            --mutations="!{osMutationsFile}" \
            --layout="!{layoutFile}" \
            --error-rates="!{onTargetErrorRatesFile}" \
            --off-target-error-rates="!{offTargetErrorRatesFile}" \
            --size-characterisation="!{sizeCharacterisationFile}" \
            --invar-scores="!{invarScoresFile}" \
            --study="!{params.STUDY}" \
            --error-suppression="!{params.ERROR_SUPPRESSION_NAME}" \
            --family-size="!{params.FAMILY_SIZE}" \
            --outlier-suppression=!{params.OUTLIER_SUPPRESSION_THRESHOLD}
        """
}

workflow analysis
{
    take:
        osMutationsChannel
        layoutChannel
        backgroundErrorRatesChannel
        offTargetErrorRatesChannel
        sizeCharacterisationChannel
        invarScoresChannel

    main:
        runAnalysis(osMutationsChannel, layoutChannel,
                    backgroundErrorRatesChannel, offTargetErrorRatesChannel,
                    sizeCharacterisationChannel, invarScoresChannel)

        generateReport(osMutationsChannel, layoutChannel,
                       backgroundErrorRatesChannel, offTargetErrorRatesChannel,
                       sizeCharacterisationChannel, invarScoresChannel)

    emit:
        plots = runAnalysis.out.plots
}
