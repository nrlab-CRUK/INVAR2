process runAnalysis
{
    memory '2g'
    time   '10m'

    publishDir params.ANALYSIS_DIR, mode: 'link', overwrite: true

    input:
        path osMutationsFile
        path tumourMutationsFile
        path layoutFile
        path onTargetErrorRatesFile
        path offTargetErrorRatesFile
        path sizeCharacterisationFile
        path invarScoresFile

    output:
        path "*.pdf", emit: "plots"
        path "${params.STUDY}_invar2_analysis.html", emit: "report"
        path "mutations_tracking.csv", emit: "mutationsTracking"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/5_analysis/analysis.R" \
            --mutations="!{osMutationsFile}" \
            --tumour-mutations="!{tumourMutationsFile}" \
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
        tumourMutationsChannel
        layoutChannel
        backgroundErrorRatesChannel
        offTargetErrorRatesChannel
        sizeCharacterisationChannel
        invarScoresChannel

    main:
        runAnalysis(osMutationsChannel, tumourMutationsChannel, layoutChannel,
                    backgroundErrorRatesChannel, offTargetErrorRatesChannel,
                    sizeCharacterisationChannel, invarScoresChannel)

    emit:
        plots = runAnalysis.out.plots
}
