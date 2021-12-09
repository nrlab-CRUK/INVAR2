process FilterMutationsAndCalculateErrorRates
{
    memory '8g'
    cpus 1
    time '1h'

    publishDir 'invar3', mode: 'link'

    input:
        path mutationFile
        path tumourMutationsFile
        path layoutFile

    output:
        path onTargetFile, emit: "onTargetFile"
        path "*.rds", emit: "rdsFiles"

    shell:
        tapasSetting = "${params.ERROR_SUPPRESSION_NAME}_BQ_${params.BASEQ}.MQ_${params.MAPQ}"
        onTargetFile = "${tapasSetting}.on_target.tsv"

        """
        export INVAR_HOME="!{projectDir}"

        Rscript --vanilla "!{projectDir}/R/invar3/invar3.R" \
            --tapas="!{tapasSetting}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            "!{mutationFile}"
        """
}


workflow invar3
{
    take:
        mutation_channel

    main:
        tumour_mutations_channel = channel.fromPath(params.TUMOUR_MUTATIONS_CSV, checkIfExists: true)
        layout_channel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

        FilterMutationsAndCalculateErrorRates(mutation_channel, tumour_mutations_channel, layout_channel)

    emit:
        onTargetFile = FilterMutationsAndCalculateErrorRates.out.onTargetFile
}
