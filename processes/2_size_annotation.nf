@Grab('com.xlson.groovycsv:groovycsv:1.3')

/**
 * Helper function for turning the output from annotateMutationsWithFragmentSize
 * into a channel containing tuples of pool, barcode, patient mutation belongs to and
 * mutations table file.
 *
 * The index file from the process needs to be read in to get the information per file.
 * That can be used to produce the list of tuples required for the downstream channel.
 * The function should be called within a closure using Nextflow's flatMap method.
 *
 * The index file has columns for the pool, barcode, patient and file name, so it's
 * a matter of finding the row that has the file name.
 */
def combineIndexFileWithMutationsFiles(indexFile, mutationsFiles)
{
    def indexContent = com.xlson.groovycsv.CsvParser.parseCsv(indexFile.getText('UTF-8'))

    def tuples = []
    for (file in mutationsFiles)
    {
        def fileInfo = indexContent.find { row -> row.FILE_NAME == file.name }
        assert fileInfo : "No information in ${indexFile.name} found for ${file.name}"
        tuples << tuple(fileInfo.POOL, fileInfo.BARCODE, fileInfo.PATIENT_MUTATION_BELONGS_TO, file)
    }
    return tuples
}

process createSNVList
{
    executor 'local'
    memory '256m'
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
    memory '1g'

    tag "${pool} ${barcode}"

    input:
        tuple val(pool), val(barcode), path(bamFile), path(bamIndex)
        each path(snvList)

    output:
        tuple val(pool), val(barcode), path(insertsFile)

    shell:
        insertsFile = "${pool}.${barcode}.inserts.tsv"

        template "2_size_annotation/getFragmentSize.sh"
}

process annotateMutationsWithFragmentSize
{
    tag "${pool} ${barcode}"

    memory '2g'
    cpus   { Math.min(params.MAX_CORES, 4) }

    input:
        tuple val(pool), val(barcode), path(fragmentSizesFile)
        each path(mutationsFile)

    output:
        tuple val(pool), val(barcode), path("index.csv"), path("mutation_table.with_sizes.*.rds"), emit: "mutationsFiles"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/2_size_annotation/sizeAnnotation.R" \
            --mutations="!{mutationsFile}" \
            --fragment-sizes="!{fragmentSizesFile}" \
            --pool="!{pool}" --barcode="!{barcode}" \
            --threads=!{task.cpus} \
            !{params.containsKey('SAMPLING_SEED') ? "--sampling-seed=${params['SAMPLING_SEED']}" : ""}
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
         * annotateMutationsWithFragmentSize, which is a tuple of pool, barcode,
         * and index file and multiple RDS files (one per patient) into a channel
         * of pool, barcode, patient and one RDS file. The index file is used
         * as a mapping between the patient identifier and the file relevant to
         * that patient. This is better (in my opinion) that trying to extract
         * that information from the file name.
         */

        perPatientChannel =
            annotateMutationsWithFragmentSize.out.mutationsFiles
                .flatMap {
                    pool, barcode, indexFile, mutationsFiles ->
                    combineIndexFileWithMutationsFiles(indexFile, mutationsFiles)
                }

    emit:
        mutationsFiles = perPatientChannel
}
