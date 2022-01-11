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
    memory '1g'
    cpus   1
    time   '4h'

    publishDir 'insert_sizes', mode: 'link'

    input:
        tuple val(pool), val(barcode), path(bamFile), path(bamIndex)
        each path(snvList)

    output:
        tuple val(pool), val(barcode), path(insertsFile)

    shell:
        insertsFile = "${bamFile.name}.inserts_for_annotation.tsv"

        template "2_size_annotation/getFragmentSize.sh"
}


workflow sizeAnnotation
{
    main:
        patient_list_channel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)

        bam_channel = channel.fromPath(params.INPUT_FILES, checkIfExists: true)
            .splitCsv(header: true, by: 1, strip: true)
            .map {
                row ->
                bam = file(row.FILE_NAME, checkIfExists: true)
                index = file("${bam.name}.bai") // The index file may or may not exist.
                tuple row.POOL, row.BARCODE, bam, index
            }

        createSNVList(patient_list_channel)

        getFragmentSize(bam_channel, createSNVList.out)
}
