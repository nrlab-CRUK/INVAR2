include { logException } from '../functions/debugging'

import groovy.json.JsonOutput


process SlopBED
{
    executor 'local'
    memory '32m'
    cpus 1
    time '5m'
    
    input:
        path bedFile
        path genomeFile
    
    output:
        path filteredBedFile, emit: "sloppedBed"
    
    shell:
        filteredBedFile = "slopped.filtered.bed"
        
        template "invar1/slopbed.sh"
}


process mpileup
{
    cpus 1
    time '4h'
    
    publishDir "mpileup", mode: 'link'
    
    input:
        path sloppedBedFile
        path fastaReference
        each path(bamFile)
    
    output:
        path vcfFile, emit: "vcfFile"
    
    shell:
        vcfFile = "${bamFile.baseName}.vcf"
        
        dedupFlags = params.REMOVE_DUPLICATES ? "-R --ff UNMAP" : ""

        template "invar1/mpileup.sh"
}


process biallelic
{
    executor 'local'
    cpus 1
    time '5m'
    
    input:
        path vcfFile
    
    output:
        tuple val(slx), val(barcode), path(mutationFile)
    
    shell:
        slx = ""
        
        def matcher = vcfFile.name =~ /(SLX-\d+)/
        if (matcher.find())
        {
            slx = matcher[0][1]
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
    
        mutationFile = "${slx}_${barcode}.BQ_${params.BASEQ}.MQ_${params.MAPQ}.final.tsv"
        
        template "invar1/biallelic.sh"
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
        tuple val(slx), val(barcode), path(mutationFile) 
    
    output:
        tuple val(slx), val(barcode), path(tabixFile)
        
    shell:
        tabixFile = "${slx}_${barcode}_snp.vcf"
        
        template "invar1/tabix.sh"
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
        tuple val(slx), val(barcode), path(mutationFile) 
    
    output:
        tuple val(slx), val(barcode), path(tabixFile)
        
    shell:
        tabixFile = "${slx}_${barcode}_cosmic.vcf"
        
        template "invar1/tabix.sh"
}

process trinucleotide
{
    memory '256m'
    cpus 1
    time '30m'
    
    input:
        each path(fastaReference)
        tuple val(slx), val(barcode), path(mutationFile) 
    
    output:
        tuple val(slx), val(barcode), path(trinucleotideFile)
        
    shell:
        trinucleotideFile = "${slx}_${barcode}_trinucleotide.fa"
        
        template "invar1/samtools_faidx.sh"
}

process annotateMutation
{
    executor 'local'
    memory '1G'
    cpus 1
    time '1h'
    
    publishDir 'mutations', mode: 'link'
    
    input:
        tuple val(slx), val(barcode), path(mutationFile), path(snp), path(cosmic), path(trinucleotide)
    
    output:
        tuple val(slx), val(barcode), path(annotationFile)
    
    shell:
        annotationFile = "${slx}_${barcode}.mutations.tsv"
        
        """
        python3 "!{projectDir}/python/invar1/addTabixAndTrinuc.py" \
            !{mutationFile} \
            !{snp} \
            !{cosmic} \
            !{trinucleotide} \
            !{annotationFile}
        """
}

workflow invar1
{
    main:
        bed_channel = channel.fromPath(params.BED, checkIfExists: true)
        genome_channel = channel.fromPath(params.HG19_GENOME, checkIfExists: true)
        fasta_channel = channel.fromPath(params.FASTAREF, checkIfExists: true)
        snp_channel = channel.fromPath(params.K1G_DB, checkIfExists: true)
        snp_index_channel = channel.fromPath("${params.K1G_DB}.tbi")
        cosmic_channel = channel.fromPath(params.COSMIC_DB, checkIfExists: true)
        cosmic_index_channel = channel.fromPath("${params.COSMIC_DB}.tbi")
        
        SlopBED(bed_channel, genome_channel)
        
        bamList = file(params.INPUT_FILES, checkIfExists: true).readLines()
        bam_channel = channel.fromList(bamList).map { f -> file(f, checkIfExists: true) }
        
        mpileup(SlopBED.out, fasta_channel, bam_channel) | biallelic
        
        mutation_channel = biallelic.out.filter { slx, bc, f -> f.countLines() > 1 }
        
        tabixSnp(snp_channel, snp_index_channel, mutation_channel)
        tabixCosmic(cosmic_channel, cosmic_index_channel, mutation_channel)
        trinucleotide(fasta_channel, mutation_channel)
        
        by_bam_mutation_channel = mutation_channel 
            .join(tabixSnp.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)
            .join(tabixCosmic.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)
            .join(trinucleotide.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)
    
        all_mutations_channel =
            annotateMutation(by_bam_mutation_channel)
                .map { slx, barcode, file -> file }
                .collectFile(name: "${params.FINAL_PREFIX}.combined.final.ann.tsv",
                             storeDir: 'mutations', keepHeader: true, skip: 1, sort: false)

    emit:
        mutationFile = all_mutations_channel
}
