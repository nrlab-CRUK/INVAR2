process generalisedLikelihoodRatioTest
{
    memory '4g'
    cpus   { Math.min(Math.ceil(params.MAX_CORES / 2.0) as int, 10) }

    input:
        tuple val(pool), val(barcode), path(osMutationsFile)
        each path(sizeCharacterisationFile)

    output:
        tuple val(pool), val(barcode), path(invarScoresFile), emit: "invarScoresFile"
        path invarScoresTSV, emit: "invarScoresTSV"

    shell:
        invarScoresFile = "invar_scores.${pool}_${barcode}.rds"
        invarScoresTSV = "invar_scores.${pool}_${barcode}.tsv"

        """
        Rscript --vanilla "!{params.projectHome}/R/4_detection/generalisedLikelihoodRatioTest.R" \
            --threads=!{task.cpus} \
            --pool="!{pool}" --barcode="!{barcode}" \
            --mutations="!{osMutationsFile}" \
            --size-characterisation="!{sizeCharacterisationFile}" \
            !{params.is_bloodspot ? "--bloodspot" : ""} \
            --outlier-suppression=!{params.outlier_suppression_threshold} \
            --minimum-fragment-length=!{params.minimum_fragment_length} \
            --maximum-fragment-length=!{params.maximum_fragment_length} \
            --smoothing=!{params.SMOOTH} \
            !{params.only_weigh_mutants ? "--only-weigh-mutants" : ""} \
            !{params.containsKey('sampling_seed') ? "--sampling-seed=${params['sampling_seed']}" : ""}
        """
}

process combineGeneralisedLikelihoodRatioTestResults
{
    memory '4g'
    time   '5m'

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
            generalisedLikelihoodRatioTest.out.invarScoresFile
                .map { p, b, f -> f }
                .combine()

        combineGeneralisedLikelihoodRatioTestResults(scoresChannel)
}
