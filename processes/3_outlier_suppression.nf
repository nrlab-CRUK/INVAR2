include { makeSafeForFileName } from '../functions/naming'

process markOutliers
{
    tag "${sampleId} ${patientMutationBelongsTo}"

    memory '6g'

    input:
        tuple val(sampleId), val(patientMutationBelongsTo), path(mutationsFile)

    output:
        tuple val(sampleId), val(patientMutationBelongsTo), path(outlierMarkedFile), emit: "mutationsFile"

    shell:
        outlierMarkedFile = "mutation_table.outliersuppressed." + makeSafeForFileName(sampleId) + '.' +
                            makeSafeForFileName(patientMutationBelongsTo) + ".rds"

        """
        Rscript --vanilla "!{params.projectHome}/R/3_outlier_suppression/outlierSuppression.R" \
            --mutations="!{mutationsFile}" \
            --sample="!{sampleId}" \
            --patient="!{patientMutationBelongsTo}" \
            --outlier-suppression=!{params.OUTLIER_SUPPRESSION_THRESHOLD} \
            --allele-frequency-threshold=!{params.ALLELE_FREQUENCY_THRESHOLD} \
            --maximum-mutant-reads=!{params.MAXIMUM_MUTANT_READS}
        """
}

process sizeCharacterisation
{
    memory '10g'
    cpus   { Math.min(Math.ceil(params.MAX_CORES / 2.0) as int, osMutationsFiles.size()) }

    publishDir params.RESULTS_DIR, mode: 'link', overwrite: true, pattern: 'size_characterisation.rds'

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
    memory '8g'
    cpus   { Math.min(Math.ceil(params.MAX_CORES / 2.0) as int, osMutationsFiles.size()) }

    publishDir params.RESULTS_DIR, mode: 'link', overwrite: true

    input:
        path mutationsFile
        path osMutationsFiles

    output:
        path 'mutation_table.rds', emit: "mutationsFile"

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
                .map { s, pt, mfile -> mfile }
                .collect()

        sizeCharacterisation(outlierMarkedFiles)

        annotateMutationsWithOutlierSuppression(mutationsChannel, outlierMarkedFiles)

    emit:
        mutationsFile = annotateMutationsWithOutlierSuppression.out.mutationsFile
        sizeCharacterisationFile = sizeCharacterisation.out.summaryFile
        perSampleMutationsFiles = markOutliers.out.mutationsFile
}
