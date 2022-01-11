process slopPatientInfo
{
    executor 'local'
    memory '256m'
    cpus 1
    time '5m'

    input:
        path csvFile
        path genomeFile

    output:
        path filteredBedFile, emit: "sloppedBed"

    shell:
        filteredBedFile = "${csvFile.baseName}.slopped.filtered.bed"

        template "1_parse/slopPatientInfo.sh"
}


process mpileup
{
    cpus 1
    time '4h'

    input:
        path sloppedBedFile
        path fastaReference
        each path(bamFile)

    output:
        path vcfFile, emit: "vcfFile"

    shell:
        vcfFile = "${bamFile.baseName}.vcf"

        dedupFlags = params.REMOVE_DUPLICATES ? "-R --ff UNMAP" : ""

        template "1_parse/mpileup.sh"
}


process biallelic
{
    executor 'local'
    cpus 1
    time '5m'

    input:
        path vcfFile

    output:
        tuple val(pool), val(barcode), path(mutationFile)

    shell:
        pool = ""

        def matcher = vcfFile.name =~ /(SLX-\d+)/
        if (matcher.find())
        {
            pool = matcher[0][1]
        }

        switch (params.LIBRARY_PREP)
        {
            case "Rubicon":
                matcher = vcfFile.name =~ /(D.*_D.*)\./
                break

            case "Aglient":
                matcher = vcfFile.name =~ /(SXT[A-Z]{2}[0-9]{3})/
                break

            default:
                throw new Exception("Library prep '${params.LIBRARY_PREP}' not recognised.")
        }

        barcode = ""
        if (matcher.find())
        {
            barcode = matcher[0][1]
        }

        mutationFile = "${pool}_${barcode}.BQ_${params.BASEQ}.MQ_${params.MAPQ}.final.tsv"

        template "1_parse/biallelic.sh"
}


process tabixSnp
{
    executor 'local'
    memory '32m'
    cpus 1
    time '1h'

    input:
        each path(tabixDatabase)
        each path(tabixDatabaseIndex)
        tuple val(pool), val(barcode), path(mutationFile)

    output:
        tuple val(pool), val(barcode), path(tabixFile)

    shell:
        tabixFile = "${pool}_${barcode}_snp.vcf"

        template "1_parse/tabix.sh"
}

process tabixCosmic
{
    executor 'local'
    memory '32m'
    cpus 1
    time '1h'

    input:
        each path(tabixDatabase)
        each path(tabixDatabaseIndex)
        tuple val(pool), val(barcode), path(mutationFile)

    output:
        tuple val(pool), val(barcode), path(tabixFile)

    shell:
        tabixFile = "${pool}_${barcode}_cosmic.vcf"

        template "1_parse/tabix.sh"
}

process trinucleotide
{
    memory '256m'
    cpus 1
    time '30m'

    input:
        each path(fastaReference)
        tuple val(pool), val(barcode), path(mutationFile)

    output:
        tuple val(pool), val(barcode), path(trinucleotideFile)

    shell:
        trinucleotideFile = "${pool}_${barcode}_trinucleotide.fa"

        template "1_parse/samtools_faidx.sh"
}

process annotateMutation
{
    executor 'local'
    memory '1G'
    cpus 1
    time '1h'

    input:
        tuple val(pool), val(barcode), path(mutationFile), path(snp), path(cosmic), path(trinucleotide)

    output:
        tuple val(pool), val(barcode), path(annotatedFile)

    shell:
        annotatedFile = "${pool}_${barcode}.mutations.tsv"

        """
        python3 "!{projectDir}/python/1_parse/addTabixAndTrinucleotides.py" \
            !{mutationFile} \
            !{snp} \
            !{cosmic} \
            !{trinucleotide} \
            !{annotatedFile}
        """
}

process combineCSV
{
    executor 'local'
    memory '32m'
    cpus 1
    time '1h'

    publishDir 'unprocessed_mutations', mode: 'link'

    input:
        path(csvFiles)

    output:
        path(combinedFile)

    shell:
        combinedFile = "${params.FINAL_PREFIX}.combined.final.ann.tsv"

        template "1_parse/catCSV.sh"
}

workflow invar12
{
    take:
        bamChannel
        tumourMutationsChannel

    main:
        simpleBamChannel = bamChannel.map { pool, barcode, bam, index -> bam }

        genome_channel = channel.fromPath(params.HG19_GENOME, checkIfExists: true)
        fasta_channel = channel.fromPath(params.FASTAREF, checkIfExists: true)
        snp_channel = channel.fromPath(params.K1G_DB, checkIfExists: true)
        snp_index_channel = channel.fromPath("${params.K1G_DB}.tbi")
        cosmic_channel = channel.fromPath(params.COSMIC_DB, checkIfExists: true)
        cosmic_index_channel = channel.fromPath("${params.COSMIC_DB}.tbi")

        slopPatientInfo(tumourMutationsChannel, genome_channel)

        mpileup(slopPatientInfo.out, fasta_channel, simpleBamChannel) | biallelic

        mutationChannel = biallelic.out.filter { pool, bc, f -> f.countLines() > 1 }

        tabixSnp(snp_channel, snp_index_channel, mutationChannel)
        tabixCosmic(cosmic_channel, cosmic_index_channel, mutationChannel)
        trinucleotide(fasta_channel, mutationChannel)

        byBamMutationChannel = mutationChannel
            .join(tabixSnp.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)
            .join(tabixCosmic.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)
            .join(trinucleotide.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)

        allMutationsChannel =
            annotateMutation(byBamMutationChannel)
                .collect { pool, barcode, file -> file }

        combineCSV(allMutationsChannel)

    emit:
        mutationFile = combineCSV.out
}
