include { makeSafeForFileName } from '../functions/naming'

process generalisedLikelihoodRatioTest
{
    // This process is given up to ten cores because for non patient specific
    // combinations it will run ten iterations of the algorithm, which can be
    // run in parallel.

    tag "${pool} ${barcode} ${patientMutationBelongsTo}"

    memory '4g'
    cpus   { Math.min(params.MAX_CORES, 10) }

    input:
        tuple val(pool), val(barcode), val(patientMutationBelongsTo), path(osMutationsFile)
        each path(sizeCharacterisationFile)

    output:
        tuple val(pool), val(barcode), val(patientMutationBelongsTo), path(invarScoresFile), emit: "invarScores"

    shell:
        invarScoresFile = "invar_scores.${pool}.${barcode}.${makeSafeForFileName(patientMutationBelongsTo)}.rds"

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
    memory '4g'
    time   '5m'

    publishDir 'invarScores', mode: 'link', overwrite: true

    input:
        path invarScoresFiles

    output:
        path 'invar_scores.rds', emit: 'invarScoresFile'
        path 'invar_scores.tsv', emit: 'invarScoresTSV'

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
                .map { p, b, pt, f -> f }
                .collect()

        combineGeneralisedLikelihoodRatioTestResults(scoresChannel)
}
