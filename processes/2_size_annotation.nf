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

    input:
        tuple val(pool), val(barcode), path(bamFile), path(bamIndex)
        each path(snvList)

    output:
        tuple val(pool), val(barcode), path(insertsFile)

    shell:
        insertsFile = "${bamFile.name}.inserts_for_annotation.tsv"

        template "2_size_annotation/getFragmentSize.sh"
}

process annotateMutationsWithFragmentSize
{
    memory '2g'
    cpus   1
    time   '4h'

    input:
        tuple val(pool), val(barcode), path(fragmentSizesFile)
        each path(mutationsFile)

    output:
        tuple val(pool), val(barcode), path(sampleSpecificMutationsFile), emit: "mutationsFile"
        path sampleSpecificMutationsTSV, emit: "mutationsTSV"

    shell:
        sampleSpecificMutationsFile = "mutation_table.with_sizes.${pool}_${barcode}.rds"
        sampleSpecificMutationsTSV = "mutation_table.with_sizes.${pool}_${barcode}.tsv"

        """
        Rscript --vanilla "!{params.projectHome}/R/2_size_annotation/sizeAnnotation.R" \
            --mutations="!{mutationsFile}" \
            --fragment-sizes="!{fragmentSizesFile}" \
            --pool="!{pool}" --barcode="!{barcode}" \
            !{params.containsKey('sampling_seed') ? "--sampling-seed=${params.sampling_seed}" : ""}
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

    emit:
        mutationsFiles = annotateMutationsWithFragmentSize.out.mutationsFile
}