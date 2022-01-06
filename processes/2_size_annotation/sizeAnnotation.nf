process createSNVList
{
    executor 'local'
    memory '256m'
    cpus 1
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
    memory '2g'
    cpus   1
    time   '4h'

    publishDir 'insert_sizes', mode: 'link'

    input:
        tuple path(bamFile), path(bamIndex)
        each path(snvList)

    output:
        path insertsFile

    shell:
        insertsFile = "${bamFile.name}.inserts_for_annotation.csv"

        template "2_size_annotation/getFragmentSize.sh"
}


workflow sizeAnnotation
{
    main:
        patient_list_channel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)

        bamList = file(params.INPUT_FILES, checkIfExists: true).readLines()
        bam_channel = channel.fromList(bamList)
            .map {
                name ->
                f = file(name, checkIfExists: true)
                index = file("${name}.bai") // The index file may or may not exist.
                tuple f, index
            }

        createSNVList(patient_list_channel)

        getFragmentSize(bam_channel, createSNVList.out)
}
