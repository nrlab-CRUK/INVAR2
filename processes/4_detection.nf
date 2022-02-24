include { makeSafeForFileName } from '../functions/naming'

process generalisedLikelihoodRatioTest
{
    // This process is given up to ten cores because for non patient specific
    // combinations it will run ten iterations of the algorithm, which can be
    // run in parallel.

    tag "${sampleId} ${patientMutationBelongsTo}"

    memory '48g'
    cpus   { Math.min(params.MAX_CORES, 10) }

    input:
        tuple val(sampleId), val(patientMutationBelongsTo), path(osMutationsFile)
        each path(sizeCharacterisationFile)

    output:
        tuple val(sampleId), val(patientMutationBelongsTo), path(invarScoresFile), emit: "invarScores"

    shell:
        invarScoresFile = "invar_scores." + makeSafeForFileName(sampleId) + '.' + makeSafeForFileName(patientMutationBelongsTo) + ".rds"

        """
        Rscript --vanilla "!{params.projectHome}/R/4_detection/generalisedLikelihoodRatioTest.R" \
            --threads=!{task.cpus} \
            --mutations="!{osMutationsFile}" \
            --size-characterisation="!{sizeCharacterisationFile}" \
            !{params.IS_BLOODSPOT ? "--bloodspot" : ""} \
            --outlier-suppression=!{params.OUTLIER_SUPPRESSION_THRESHOLD} \
            --minimum-fragment-length=!{params.MINIMUM_FRAGMENT_LENGTH} \
            --maximum-fragment-length=!{params.MAXIMUM_FRAGMENT_LENGTH} \
            --smoothing=!{params.SMOOTHING} \
            !{params.ONLY_WEIGH_MUTANTS ? "--only-weigh-mutants" : ""} \
            !{params.containsKey('SAMPLING_SEED') ? "--sampling-seed=${params['SAMPLING_SEED']}" : ""}
        """
}

process combineGeneralisedLikelihoodRatioTestResults
{
    memory '256m'
    time   '5m'

    publishDir params.RESULTS_DIR, mode: 'link', overwrite: true

    input:
        path invarScoresFiles

    output:
        path 'invar_scores.rds', emit: 'invarScores'

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/4_detection/combineGeneralisedLikelihoodRatioTestResults.R" \
            !{invarScoresFiles}
        """

}


workflow detection
{
    take:
        osMutationsChannel
        sizeCharacterisationChannel

    main:
        generalisedLikelihoodRatioTest(osMutationsChannel, sizeCharacterisationChannel)

        scoresChannel =
            generalisedLikelihoodRatioTest.out.invarScores
                .map { s, pt, f -> f }
                .collect()

        combineGeneralisedLikelihoodRatioTestResults(scoresChannel)

    emit:
        invarScores = combineGeneralisedLikelihoodRatioTestResults.out.invarScores
}
