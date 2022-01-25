process markOutliers
{
    memory '2g'
    cpus   1
    time   '4h'

    input:
        tuple val(pool), val(barcode), path(mutationsFile)

    output:
        tuple val(pool), val(barcode), path(outlierMarkedFile), emit: "mutationsFile"
        path outlierMarkedTSV, emit: "mutationsTSV"

    shell:
        outlierMarkedFile = "mutation_table.outliersuppressed.${pool}_${barcode}.rds"
        outlierMarkedTSV = "mutation_table.outliersuppressed.${pool}_${barcode}.tsv"

        """
        Rscript --vanilla "!{params.projectHome}/R/3_outlier_suppression/outlierSuppression.R" \
            --mutations="!{mutationsFile}" \
            --pool="!{pool}" --barcode="!{barcode}" \
            --outlier-suppression=!{params.outlier_suppression_threshold}
        """
}

process sizeCharacterisation
{
    memory '4g'
    cpus   6
    time   '1h'

    input:
        path mutationsFiles

    output:
        path 'size_characterisation.all.rds', emit: "allSizesFile"
        path 'size_characterisation.all.tsv', emit: "allSizesTSV"
        path 'size_characterisation.rds', emit: "summaryFile"
        path 'size_characterisation.tsv', emit: "summaryTSV"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/3_outlier_suppression/sizeCharacterisation.R" \
            --threads=!{task.cpus} \
            !{mutationsFiles}
        """
}


workflow outlierSuppression
{
    take:
        mutationsChannel

    main:
        markOutliers(mutationsChannel)

        sizedFiles = markOutliers.out.mutationsFile.map { pool, barcode, mfile -> mfile }.collect()

        sizeCharacterisation(sizedFiles)

    emit:
        mutationsFiles = markOutliers.out.mutationsFile
}
