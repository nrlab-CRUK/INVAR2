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
        
        dedupFlags = params.removeDuplicates ? "-R --ff UNMAP" : ""

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
        path mutationFile, emit: "mutationFile"
    
    shell:
        mutationFile = "${vcfFile.baseName}.BQ_${params.BASEQ}.MQ_${params.MAPQ}.final.tsv"
        
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
    
        template "invar1/biallelic.sh"
}


process tabixSnp
{
    executor 'local'
    memory '32m'
    cpus 1
    time '5m'
    
    input:
        each path(snpDatabase)
        each path(snpDatabaseIndex)
        tuple val(slx), val(barcode), val(chromosome), val(position), val(mutation) 
    
    output:
        tuple val(slx), val(barcode), val(chromosome), val(position), val(mutation), path('snp.vcf')
        
    shell:
        chromosome_short = mutation['CHROM'].toUpperCase().replaceAll("CHR", "")
        position = mutation['POS'] as long
        
        mutationPosition = "${chromosome_short}:${position}-${position}"
        
        """
        tabix "!{snpDatabase}" "!{mutationPosition}" > snp.vcf
        """
}

process tabixCosmic
{
    executor 'local'
    memory '32m'
    cpus 1
    time '5m'
    
    input:
        each path(cosmicDatabase)
        each path(cosmicDatabaseIndex)
        tuple val(slx), val(barcode), val(chromosome), val(position), val(mutation) 
    
    output:
        tuple val(slx), val(barcode), val(chromosome), val(position), path('cosmic.vcf')
        
    shell:
        chromosome_short = mutation['CHROM'].toUpperCase().replaceAll("CHR", "")
        position = mutation['POS'] as long
        
        mutationPosition = "${chromosome_short}:${position}-${position}"
        
        """
        tabix "!{cosmicDatabase}" "!{mutationPosition}" > cosmic.vcf
        """
}

process trinucleotide
{
    memory '256m'
    cpus 1
    time '30m'
    
    input:
        each path(fastaReference)
        tuple val(slx), val(barcode), val(chromosome), val(position), val(mutation) 
    
    output:
        tuple val(slx), val(barcode), val(chromosome), val(position), path('trinucleotide.fa')
        
    shell:
        chromosome_short = mutation['CHROM'].toUpperCase().replaceAll("CHR", "")
        position = mutation['POS'] as long
        
        trinucleotidePosition = "chr${chromosome_short}:${position - 1}-${position + 1}"
        
        """
        samtools faidx "!{fastaReference}" "!{trinucleotidePosition}" > trinucleotide.fa
        """
}

process annotateMutation
{
    executor 'local'
    memory '32m'
    cpus 1
    time '5m'
    
    input:
        tuple val(slx), val(barcode), val(chromosome), val(position), val(mutation), path(snp), path(cosmic), path(trinucleotide)
    
    output:
        tuple val(slx), val(barcode), val(chromosome), val(position), path(annotationFile)
    
    shell:
        mutationJson = JsonOutput.toJson(mutation)
        annotationFile = "${slx}_${barcode}.mutations.tsv"
        
        """
        python3 "!{projectDir}/python/invar1/parseTabix.py" \
            '!{mutationJson}' \
            !{snp} \
            !{cosmic} \
            !{trinucleotide} \
            !{annotationFile}
        """
}

workflow invar1
{
    bed_channel = channel.fromPath(params.BED)
    genome_channel = channel.fromPath(params.HG19_GENOME)
    fasta_channel = channel.fromPath(params.FASTAREF)
    snp_channel = channel.fromPath(params.K1G_DB)
    snp_index_channel = channel.fromPath("${params.K1G_DB}.tbi")
    cosmic_channel = channel.fromPath(params.COSMIC_DB)
    cosmic_index_channel = channel.fromPath("${params.COSMIC_DB}.tbi")
    
    SlopBED(bed_channel, genome_channel)
    
    bam_channel = channel.fromPath("${launchDir}/bam/*.bam")
    
    mpileup(SlopBED.out, fasta_channel, bam_channel) | biallelic
    
    mutation_channel = biallelic.out
        .filter(f -> f.countLines() > 1)
        .splitCsv(header: true, by: 1, sep: '\t')
        .map(m -> tuple(m['SLX'], m['BARCODE'], m['CHROM'], m['POS'], m))
    
    tabixSnp(snp_channel, snp_index_channel, mutation_channel)
    tabixCosmic(cosmic_channel, cosmic_index_channel, mutation_channel)
    trinucleotide(fasta_channel, mutation_channel)
    
    combined_channel = tabixSnp.out
        .join(tabixCosmic.out, by: 0..3, failOnDuplicate: true, failOnMismatch: true)
        .join(trinucleotide.out, by: 0..3, failOnDuplicate: true, failOnMismatch: true)

    annotateMutation(combined_channel)
        .map { slx, barcode, chr, pos, file -> file }
        .collectFile(storeDir: 'mutations', keepHeader: true, skip: 1, newLine: true)
}