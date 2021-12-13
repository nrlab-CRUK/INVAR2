process createMutationsTable
{
    memory '8g'
    cpus 1
    time '1h'

    publishDir 'invar34', mode: 'link'

    input:
        path mutationFile
        path tumourMutationsFile
        path layoutFile

    output:
        path 'mutation_table.all.rds', emit: "allMutationFile"
        path 'mutation_table.filtered.rds', emit: "filteredMutationFile"

    shell:
        tapasSetting = "${params.ERROR_SUPPRESSION_NAME}_BQ_${params.BASEQ}.MQ_${params.MAPQ}"

        """
        Rscript --vanilla "!{projectDir}/R/invar34/createMutationsTable.R" \
            --mutations="!{mutationFile}" \
            --tapas="!{tapasSetting}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            --cosmic-threshold=!{params.cosmic_threshold} \
            --mqsb-threshold=!{params.individual_MQSB_threshold} \
            --max-dp=!{params.max_DP} \
            --min-ref-dp=!{params.min_ref_DP} \
            --alt-allele-threshold=!{params.n_alt_alleles_threshold} \
            --minor-alt-allele-threshold=!{params.minor_alt_allele_threshold}
        """
}

process offTargetErrorRate
{
    memory '8g'
    cpus 1
    time '1h'

    publishDir 'invar34', mode: 'link'

    input:
        each path(mutationRDSFile)
        each path(layoutFile)
        val cosmic

    output:
        path '*.rds'

    shell:
        """
        export INVAR_HOME="!{projectDir}"

        Rscript --vanilla "!{projectDir}/R/invar34/offTargetErrorRate.R" \
            --mutations="!{mutationRDSFile}" \
            --layout="!{layoutFile}" \
            !{cosmic ? "--cosmic" : ""}
        """
}


workflow invar3
{
    take:
        mutation_channel

    main:
        tumour_mutations_channel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
        layout_channel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

        createMutationsTable(mutation_channel, tumour_mutations_channel, layout_channel)

        offTargetErrorRate(createMutationsTable.out.filteredMutationFile, layout_channel, channel.of(true, false))

}
