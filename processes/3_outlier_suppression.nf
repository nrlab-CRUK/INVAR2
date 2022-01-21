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


workflow outlierSuppression
{
    take:
        mutationsChannel

    main:
        markOutliers(mutationsChannel)

    emit:
        mutationsFiles = markOutliers.out.mutationsFile
}
