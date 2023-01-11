process runAnalysis
{
    memory '2g'
    time   '3h'

    publishDir params.ANALYSIS_DIR, mode: 'link', overwrite: true

    input:
        path osMutationsFile
        path tumourMutationsFile
        path layoutFile
        path backgroundErrorRatesFile
        path onTargetLocusErrorRatesFile
        path offTargetErrorRatesFile
        path sizeCharacterisationFile
        path invarScoresFile

    output:
        path "*.pdf", emit: "plots"
        path "${params.STUDY}_invar2_analysis.html", emit: "report"
        path "mutations_tracking.csv", emit: "mutationsTracking"
        path "*.csv", emit: "dataFiles"

    shell:
        tapasSetting = "${params.ERROR_SUPPRESSION_NAME}_BQ_${params.BASE_QUALITY}.MQ_${params.MAPPING_QUALITY}"

        """
        Rscript --vanilla "!{params.projectHome}/R/5_analysis/analysis.R" \
            --mutations="!{osMutationsFile}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            --error-rates="!{backgroundErrorRatesFile}" \
            --on-target-locus-error-rates="!{onTargetLocusErrorRatesFile}" \
            --off-target-error-rates="!{offTargetErrorRatesFile}" \
            --size-characterisation="!{sizeCharacterisationFile}" \
            --invar-scores="!{invarScoresFile}" \
            --study="!{params.STUDY}" \
            --tapas="!{tapasSetting}" \
            --error-suppression="!{params.ERROR_SUPPRESSION_NAME}" \
            --family-size="!{params.FAMILY_SIZE}" \
            --outlier-suppression=!{params.OUTLIER_SUPPRESSION_THRESHOLD} \
            --max-background-allele-frequency=!{params.MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY} \
            --allele-frequency-threshold=!{params.ALLELE_FREQUENCY_THRESHOLD} \
            --max-mutant-reads=!{params.MAXIMUM_MUTANT_READS} \
            --min-informative-reads=!{params.MINIMUM_INFORMATIVE_READS} \
            --score-specificity=!{params.SCORE_SPECIFICITY}
            --LNP_threshold=!{params.PROPORTION_OF_CONTROLS}
        """
}


workflow analysis
{
    take:
        osMutationsChannel
        tumourMutationsChannel
        layoutChannel
        backgroundErrorRatesChannel
        onTargetLocusErrorRatesChannel
        offTargetErrorRatesChannel
        sizeCharacterisationChannel
        invarScoresChannel

    main:
        runAnalysis(osMutationsChannel, tumourMutationsChannel, layoutChannel,
                    backgroundErrorRatesChannel, onTargetLocusErrorRatesChannel,
                    offTargetErrorRatesChannel,
                    sizeCharacterisationChannel, invarScoresChannel)
}
