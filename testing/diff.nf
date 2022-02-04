process rdsDiff
{
    tag "${name} ${rdsFile.name}"

    errorStrategy 'ignore'

    input:
        tuple val(name), path(rdsFile), path(reference)
    
    shell:
        tsvFile = rdsFile.name.replaceAll(/.rds$/, /.tsv/)
        
        """
        Rscript --vanilla "!{params.projectHome}/testing/R/rdsToTSV.R" \
            "!{rdsFile}" "!{tsvFile}"

        diff "!{tsvFile}" "!{reference}"
        """
}

process diff
{
    tag "${name} ${generated.name}"

    memory '64m'

    errorStrategy 'ignore'

    input:
        tuple val(name), path(generated), path(reference)

    shell:
        """
        diff "!{generated}" "!{reference}"
        """
}
