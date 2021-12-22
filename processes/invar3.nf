process createMutationsTable
{
    memory '8g'
    cpus 1
    time '1h'

    publishDir 'mutations', mode: 'link'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile

    output:
        path 'mutation_table.all.rds', emit: "allMutationsFile", optional: true
        path 'mutation_table.filtered.rds', emit: "filteredMutationsFile"

    shell:
        tapasSetting = "${params.ERROR_SUPPRESSION_NAME}_BQ_${params.BASEQ}.MQ_${params.MAPQ}"

        """
        Rscript --vanilla "!{projectDir}/R/invar34/createMutationsTable.R" \
            --mutations="!{mutationsFile}" \
            --tapas="!{tapasSetting}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            --cosmic-threshold=!{params.cosmic_threshold} \
            --mqsb-threshold=!{params.individual_MQSB_threshold} \
            --max-dp=!{params.max_DP} \
            --min-ref-dp=!{params.min_ref_DP} \
            --alt-alleles-threshold=!{params.n_alt_alleles_threshold} \
            --minor-alt-allele-threshold=!{params.minor_alt_allele_threshold}
        """
}

process offTargetErrorRates
{
    memory '2g'
    cpus 2
    time '1h'

    publishDir 'off_target', mode: 'link'

    input:
        path mutationsFile
        path layoutFile

    output:
        path 'locus_error_rates.off_target.rds', emit: 'locusErrorRates'
        path 'mutation_table.error_rates.cosmic.rds', emit: "cosmicErrorRates"
        path 'mutation_table.error_rates.no_cosmic.rds', emit: "noCosmicErrorRates"
        path '*.tsv', optional: true

    shell:
        """
        Rscript --vanilla "!{projectDir}/R/invar34/offTargetErrorRate.R" \
            --mutations="!{mutationsFile}" \
            --layout="!{layoutFile}" \
            --control-proportion=!{params.proportion_of_controls} \
            --max-background-af=!{params.max_background_mean_AF} \
            !{params.is_bloodspot ? "--bloodspot" : ""}
        """
}

process createOnTargetMutationsTable
{
    memory '2g'
    cpus 2
    time '1h'

    publishDir 'on_target', mode: 'link'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile
        path errorRatesFile

    output:
        path 'mutation_table.on_target.rds', emit: "onTargetMutationsFile"
        path '*.tsv', optional: true

    shell:
        """
        Rscript --vanilla "!{projectDir}/R/invar34/createOnTargetMutationsTable.R" \
            --mutations="!{mutationsFile}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            --error-rates="!{errorRatesFile}"
        """
}


workflow invar3
{
    take:
        mutationsChannel

    main:
        tumourMutationsChannel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
        layoutChannel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

        createMutationsTable(mutationsChannel, tumourMutationsChannel, layoutChannel)

        offTargetErrorRates(createMutationsTable.out.filteredMutationsFile, layoutChannel)

        createOnTargetMutationsTable(
            createMutationsTable.out.filteredMutationsFile,
            tumourMutationsChannel,
            layoutChannel,
            offTargetErrorRates.out.noCosmicErrorRates)
}
