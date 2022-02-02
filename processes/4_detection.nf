@Grab('com.xlson.groovycsv:groovycsv:1.3')

/**
 * Helper function for turning the output from splitForGeneralisedLikelihoodRatioTest
 * into a channel containing tuples of pool barcode, patient mutation belongs to and
 * mutations table file.
 *
 * The index file from the process needs to be read in to get the information per file.
 * That can be used to produce the list of tuples required for the downstream channel.
 * The function should be called within a closure using Nextflow's flatMap method.
 */
def combineIndexFileWithMutationsFiles(indexFile, mutationsFiles)
{
    def indexContent = com.xlson.groovycsv.CsvParser.parseCsv(indexFile.getText('UTF-8'))
    def filenameToInfo = [:]
    for (row in indexContent)
    {
        filenameToInfo[row.FILE_NAME] = row
    }
    def tuples = []
    for (file in mutationsFiles)
    {
        def row = filenameToInfo[file.name]
        tuples << tuple(row.POOL, row.BARCODE, row.PATIENT_MUTATION_BELONGS_TO, file)
    }
    return tuples
}

/**
 * Convert a patient id into a safe version for a file name. Spaces are replaced
 * by a single underscore then any character that is not a word character (letter, digit
 * or underscore) is removed.
 */
def safePatientId(id)
{
    id.replaceAll(/\s+/, '_').replace(/[^\w]+/, '')
}


process splitForGeneralisedLikelihoodRatioTest
{
    cpus   { Math.min(Math.ceil(params.MAX_CORES / 2.0) as int, 4) }

    input:
        tuple val(pool), val(barcode), path(osMutationsFile)

    output:
        tuple val(pool), val(barcode), path("index.csv"), path("*.perpatient.rds"), emit: "perPatientFiles"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/4_detection/splitForGeneralisedLikelihoodRatioTest.R" \
            --threads=!{task.cpus} \
            --pool="!{pool}" --barcode="!{barcode}" \
            --mutations="!{osMutationsFile}"
        """
}

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
        tuple val(pool), val(barcode), val(patientMutationBelongsTo), path(invarScoresFile), emit: "invarScoresFile"
        path invarScoresTSV, emit: "invarScoresTSV"

    shell:
        invarScoresFile = "invar_scores.${pool}.${barcode}.${safePatientId(patientMutationBelongsTo)}.rds"
        invarScoresTSV = "invar_scores.${pool}.${barcode}.${safePatientId(patientMutationBelongsTo)}.tsv"

        """
        Rscript --vanilla "!{params.projectHome}/R/4_detection/generalisedLikelihoodRatioTest.R" \
            --threads=!{task.cpus} \
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

    publishDir 'invarScores'

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
        splitForGeneralisedLikelihoodRatioTest(osMutationsChannel)

        perPatientChannel =
            splitForGeneralisedLikelihoodRatioTest.out.perPatientFiles
                .flatMap {
                    pool, barcode, indexFile, mutationsFiles ->
                    combineIndexFileWithMutationsFiles(indexFile, mutationsFiles)
                }

        generalisedLikelihoodRatioTest(perPatientChannel, sizeCharacterisationChannel)

        scoresChannel =
            generalisedLikelihoodRatioTest.out.invarScoresFile
                .map { p, b, pt, f -> f }
                .collect()

        combineGeneralisedLikelihoodRatioTestResults(scoresChannel)
}
