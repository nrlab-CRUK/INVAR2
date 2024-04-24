@Grab('com.xlson.groovycsv:groovycsv:1.3')

/**
 * Helper function for turning the output from annotateMutationsWithFragmentSize
 * into a channel containing tuples of sample id, patient mutation belongs to and
 * mutations table file.
 *
 * The index file from the process needs to be read in to get the information per file.
 * That can be used to produce the list of tuples required for the downstream channel.
 * The function should be called within a closure using Nextflow's flatMap method.
 *
 * The index file has columns for the sample id, patient and file name, so it's
 * a matter of finding the row that has the file name.
 */
def combineIndexFileWithMutationsFiles(indexFile, mutationsFiles)
{
    def tuples = []
    indexFile.withReader
    {
        reader ->

        for (row in com.xlson.groovycsv.CsvParser.parseCsv(reader))
        {
            def file = mutationsFiles.find { f -> f.name == row.FILE_NAME }
            assert file : "No file found for file name ${row.FILE_NAME} listed in ${indexFile.name}"
            tuples << tuple(row.SAMPLE_ID, row.PATIENT_MUTATION_BELONGS_TO, file)
        }
    }
    return tuples
}

include { makeSafeForFileName } from '../functions/naming'


process createSNVList
{
    executor 'local'
    memory '128m'
    time '5m'

    input:
        path csvFile

    output:
        path snvListFile, emit: "snvList"

    shell:
        snvListFile = "snvlist.txt"

        template "2_size_annotation/createSNVList.sh"
}

process getFragmentSize
{
    time '24h'
    memory '30g'

    tag "${sampleId}"

    input:
        tuple val(sampleId), path(bamFile), path(bamIndex)
        each path(snvList)

    output:
        tuple val(sampleId), path(insertsFile)

    shell:
        insertsFile = makeSafeForFileName(sampleId) + ".inserts.tsv"

        template "2_size_annotation/getFragmentSize.sh"
}

process annotateMutationsWithFragmentSize
{
    tag "${sampleId}"

    memory '400g'
    cpus   { Math.min(params.MAX_CORES, 8) }

    input:
        tuple val(sampleId), path(fragmentSizesFile)
        each path(mutationsFile)

    output:
        tuple val(sampleId), path("index.csv"), path("mutation_table.with_sizes.*.rds"), emit: "mutationsFiles"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/2_size_annotation/sizeAnnotation.R" \
            --mutations="!{mutationsFile}" \
            --fragment-sizes="!{fragmentSizesFile}" \
            --sample="!{sampleId}" \
            --threads=!{task.cpus}
        """
}


workflow sizeAnnotation
{
    take:
        bamChannel
        mutationsChannel
        tumourMutationsChannel

    main:
        createSNVList(tumourMutationsChannel)

        getFragmentSize(bamChannel, createSNVList.out)

        annotateMutationsWithFragmentSize(getFragmentSize.out, mutationsChannel)

        /*
         * This bit of Nextflow jiggery pokery converts the output from
         * annotateMutationsWithFragmentSize, which is a tuple of sampleId,
         * index file and multiple RDS files (one per patient) into a channel
         * of sampleId, patient and one RDS file. The index file is used
         * as a mapping between the patient identifier and the file relevant to
         * that patient. This is better (in my opinion) that trying to extract
         * that information from the file name.
         */

        perPatientChannel =
            annotateMutationsWithFragmentSize.out.mutationsFiles
                .flatMap {
                    sampleId, indexFile, mutationsFiles ->
                    combineIndexFileWithMutationsFiles(indexFile, mutationsFiles)
                }

    emit:
        mutationsFiles = perPatientChannel
}
