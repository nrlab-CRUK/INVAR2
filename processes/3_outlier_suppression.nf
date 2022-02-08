include { makeSafeForFileName } from '../functions/naming'

process markOutliers
{
    tag "${pool} ${barcode} ${patientMutationBelongsTo}"

    memory '2g'

    input:
        tuple val(pool), val(barcode), val(patientMutationBelongsTo), path(mutationsFile)

    output:
        tuple val(pool), val(barcode), val(patientMutationBelongsTo), path(outlierMarkedFile), emit: "mutationsFile"

    shell:
        outlierMarkedFile = "mutation_table.outliersuppressed.${pool}.${barcode}.${makeSafeForFileName(patientMutationBelongsTo)}.rds"

        """
        Rscript --vanilla "!{params.projectHome}/R/3_outlier_suppression/outlierSuppression.R" \
            --mutations="!{mutationsFile}" \
            --outlier-suppression=!{params.OUTLIER_SUPPRESSION_THRESHOLD}
        """
}

process sizeCharacterisation
{
    memory '4g'
    cpus   { Math.min(Math.ceil(params.MAX_CORES / 2.0) as int, osMutationsFiles.size()) }

    input:
        path osMutationsFiles

    output:
        path 'size_characterisation.all.rds', emit: "allSizesFile"
        path 'size_characterisation.rds', emit: "summaryFile"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/3_outlier_suppression/sizeCharacterisation.R" \
            --threads=!{task.cpus} \
            !{osMutationsFiles}
        """
}

process annotateMutationsWithOutlierSuppression
{
    memory '4g'
    cpus   { Math.min(Math.ceil(params.MAX_CORES / 2.0) as int, osMutationsFiles.size()) }

    input:
        path mutationsFile
        path osMutationsFiles

    output:
        path 'mutation_table.with_outliers.rds', emit: "mutationsFile"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/3_outlier_suppression/outlierSuppressionAnnotation.R" \
            --threads=!{task.cpus} \
            --mutations="!{mutationsFile}" \
            !{osMutationsFiles}
        """
}


workflow outlierSuppression
{
    take:
        mutationsChannel
        sizedMutationsChannel

    main:
        markOutliers(sizedMutationsChannel)

        outlierMarkedFiles =
            markOutliers.out.mutationsFile
                .map { p, b, pt, mfile -> mfile }
                .collect()

        sizeCharacterisation(outlierMarkedFiles)

        annotateMutationsWithOutlierSuppression(mutationsChannel, outlierMarkedFiles)

    emit:
        mutationsFile = annotateMutationsWithOutlierSuppression.out.mutationsFile
        sizeCharacterisationFile = sizeCharacterisation.out.summaryFile
        perSampleMutationsFiles = markOutliers.out.mutationsFile
}
