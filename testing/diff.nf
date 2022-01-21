process diff
{
    tag "${name} ${generated.name}"

    memory '64m'
    cpus 1
    time '10m'

    errorStrategy 'ignore'

    input:
        tuple val(name), path(generated), path(reference)

    shell:
        """
            diff "!{generated}" "!{reference}"
        """
}
