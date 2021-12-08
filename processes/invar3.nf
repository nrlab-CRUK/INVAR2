process FilterMutationsAndCalculateRates
{
    memory '8g'
    cpus 1
    time '1h'
    
    publishDir 'invar3', mode: 'link'

    input:
        path mutationFile
        path patientBedFile
        path layoutFile

    output:
        path onTargetFile, emit: "onTargetFile"
        path "*.rds", emit: "rdsFiles"

    shell:
        tapasSetting = "${params.ERROR_SUPPRESSION_NAME}_BQ_${params.BASEQ}.MQ_${params.MAPQ}"
        onTargetFile = "${tapasSetting}.on_target.tsv"
        
        """
        Rscript --vanilla "!{projectDir}/R/invar3/invar3.R" \
            --mutations-list="!{patientBedFile}" \
            --tapas="!{tapasSetting}" \
            --layout="!{layoutFile}" \
            "!{mutationFile}"
        """
}


workflow invar3
{
    take:
        mutation_channel

    main:
        bed_channel = channel.fromPath(params.BED, checkIfExists: true)
        layout_channel = channel.fromPath(params.LAYOUT_TABLE, checkIfExists: true)

        FilterMutationsAndCalculateRates(mutation_channel, bed_channel, layout_channel)

    emit:
        onTargetFile = FilterMutationsAndCalculateRates.out.onTargetFile
}
