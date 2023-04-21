include { makeSafeForFileName } from '../functions/naming'

process createSequenceDictionary
{
    time '5m'

    input:
        path fastaReference
        path fastaIndex

    output:
        path sequenceDictionary, emit: "dictionary"

    shell:
        javaMem = task.memory.toMega() - 128
        sequenceDictionary = "${fastaReference.baseName}.dict"

        template "1_parse/createSequenceDictionary.sh"
}

process slopPatientInfo
{
    executor 'local'
    memory '128m'
    time '5m'

    input:
        path csvFile
        path dictionaryFile

    output:
        path filteredBedFile, emit: "sloppedBed"

    shell:
        filteredBedFile = "${csvFile.baseName}.slopped.filtered.bed"

        template "1_parse/slopPatientInfo.sh"
}


process mpileup
{
    tag "${sampleId}"

    cpus { Math.min(params.MAX_CORES, 2) }
    time "12h"
    

    input:
        each path(sloppedBedFile)
        each path(fastaReference)
        each path(fastaIndex)
        tuple val(sampleId), path(bamFile), path(bamIndex)

    output:
        tuple val(sampleId), path(vcfFile), emit: "vcfFile"

    shell:
        vcfFile = makeSafeForFileName(sampleId) + ".mpileup.vcf"

        dedupFlags = params.REMOVE_DUPLICATES ? "--ignore-RG --excl-flags UNMAP" : ""

        template "1_parse/mpileup.sh"
}


process biallelic
{
    tag "${sampleId}"

    memory '512m'
    time   '30m'

    input:
        tuple val(sampleId), path(vcfFile)

    output:
        tuple val(sampleId), path(mutationsFile), emit: "mutationsFile"

    shell:
        mutationsFile = makeSafeForFileName(sampleId) + ".tsv"

        template "1_parse/biallelic.sh"
}


process tabixSnp
{
    tag "${sampleId}"

    memory '64m'

    input:
        each path(tabixDatabase)
        each path(tabixDatabaseIndex)
        tuple val(sampleId), path(mutationFile)

    output:
        tuple val(sampleId), path(tabixFile)

    shell:
        tabixFile = makeSafeForFileName(sampleId) + ".snp.vcf"

        template "1_parse/tabix.sh"
}

process tabixCosmic
{
    tag "${sampleId}"

    memory '64m'

    input:
        each path(tabixDatabase)
        each path(tabixDatabaseIndex)
        tuple val(sampleId), path(mutationFile)

    output:
        tuple val(sampleId), path(tabixFile)

    shell:
        tabixFile = makeSafeForFileName(sampleId) + ".cosmic.vcf"

        template "1_parse/tabix.sh"
}

process trinucleotide
{
    tag "${sampleId}"

    memory '32m'

    input:
        each path(fastaReference)
        each path(fastaIndex)
        tuple val(sampleId), path(mutationFile)

    output:
        tuple val(sampleId), path(trinucleotideFile)

    shell:
        trinucleotideFile = makeSafeForFileName(sampleId) + ".trinucleotide.fa"

        template "1_parse/samtools_faidx.sh"
}

process annotateMutation
{
    tag "${sampleId}"

    input:
        tuple val(sampleId), path(mutationFile), path(snp), path(cosmic), path(trinucleotide)

    output:
        tuple val(sampleId), path(annotatedFile), emit: "mutationsFile"

    shell:
        annotatedFile = makeSafeForFileName(sampleId) + ".mutations.tsv"

        """
        python3 "!{params.projectHome}/python/1_parse/addTabixAndTrinucleotides.py" \
            !{mutationFile} \
            !{snp} \
            !{cosmic} \
            !{trinucleotide} \
            !{annotatedFile}
        """
}

process createMutationsTable
{
    memory '40g'
    cpus   { Math.min(params.MAX_CORES, mutationsFiles.size()) }

    input:
        path mutationsFiles
        path tumourMutationsFile

    output:
        path 'mutation_table.filtered.rds', emit: "filteredMutationsFile"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/1_parse/createMutationsTable.R" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --cosmic-threshold=!{params.COSMIC_THRESHOLD} \
            --mqsb-threshold=!{params.MQSB_THRESHOLD} \
            --max-depth=!{params.MAXIMUM_DEPTH} \
            --min-ref-depth=!{params.MINIMUM_REFERENCE_DEPTH} \
            --alt-alleles-threshold=!{params.ALT_ALLELES_THRESHOLD} \
            --minor-alt-allele-threshold=!{params.MINOR_ALT_ALLELE_THRESHOLD} \
            --threads=!{task.cpus} \
            !{mutationsFiles}
        """
}

process offTargetErrorRates
{
    memory '64g'
    cpus   { Math.min(params.MAX_CORES, 2) }

    publishDir params.RESULTS_DIR, mode: 'link', overwrite: true, pattern: "error_rates.*.rds"

    input:
        path mutationsFile
        path layoutFile

    output:
        path 'locus_error_rates.off_target.rds', emit: 'locusErrorRates'
        path 'error_rates.off_target.cosmic.rds', emit: "cosmicErrorRates"
        path 'error_rates.off_target.no_cosmic.rds', emit: "noCosmicErrorRates"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/1_parse/offTargetErrorRates.R" \
            --mutations="!{mutationsFile}" \
            --layout="!{layoutFile}" \
            --control-proportion=!{params.PROPORTION_OF_CONTROLS} \
            --max-background-allele-frequency=!{params.MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY} \
            !{params.IS_BLOODSPOT ? "--bloodspot" : ""}
        """
}

process createOnTargetMutationsTable
{
    memory '8g'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile
        path errorRatesFile

    output:
        path 'mutation_table.on_target.all.rds', emit: "onTargetMutationsFile"
        path 'background_error_rates.rds', emit: "backgroundErrorRates"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/1_parse/createOnTargetMutationsTable.R" \
            --mutations="!{mutationsFile}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            --error-rates="!{errorRatesFile}"
        """
}

process onTargetErrorRatesAndFilter
{
    memory '8g'

    publishDir params.RESULTS_DIR, mode: 'link', overwrite: true, pattern: "locus_error_rates*"

    input:
        path mutationsFile
        path tumourMutationsFile

    output:
        path 'locus_error_rates.on_target.rds', emit: 'locusErrorRates'
        path 'locus_error_rates.on_target.pdf', optional: true, emit: 'locusErrorRatesPlot'
        path 'mutation_table.on_target.rds', emit: "onTargetMutationsFile"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/1_parse/onTargetErrorRatesAndFilter.R" \
            --mutations="!{mutationsFile}" \
            --study="!{params.STUDY}" \
            --cosmic-threshold=!{params.COSMIC_THRESHOLD} \
            --control-proportion=!{params.PROPORTION_OF_CONTROLS} \
            --max-background-allele-frequency=!{params.MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY} \
            !{params.IS_BLOODSPOT ? "--bloodspot" : ""} \
            --allele-frequency-threshold=!{params.ALLELE_FREQUENCY_THRESHOLD}
        """
}


workflow parse
{
    take:
        bamChannel

        tumourMutationsChannel
        layoutChannel

    main:
        fasta_channel = channel.fromPath(params.FASTA_REFERENCE, checkIfExists: true)
        fasta_index_channel = channel.fromPath("${params.FASTA_REFERENCE}.fai")
        snp_channel = channel.fromPath(params.THOUSAND_GENOMES_DATABASE, checkIfExists: true)
        snp_index_channel = channel.fromPath("${params.THOUSAND_GENOMES_DATABASE}.tbi")
        cosmic_channel = channel.fromPath(params.COSMIC_DATABASE, checkIfExists: true)
        cosmic_index_channel = channel.fromPath("${params.COSMIC_DATABASE}.tbi")

        createSequenceDictionary(fasta_channel, fasta_index_channel)

        slopPatientInfo(tumourMutationsChannel, createSequenceDictionary.out)

        mpileup(slopPatientInfo.out, fasta_channel, fasta_index_channel, bamChannel)
        biallelic(mpileup.out)

        mutationChannel = biallelic.out.filter { sampleId, f -> f.countLines() > 1 }

        tabixSnp(snp_channel, snp_index_channel, mutationChannel)
        tabixCosmic(cosmic_channel, cosmic_index_channel, mutationChannel)
        trinucleotide(fasta_channel, fasta_index_channel, mutationChannel)

        byBamMutationChannel = mutationChannel
            .join(tabixSnp.out, by: 0, failOnDuplicate: true, failOnMismatch: true)
            .join(tabixCosmic.out, by: 0, failOnDuplicate: true, failOnMismatch: true)
            .join(trinucleotide.out, by: 0, failOnDuplicate: true, failOnMismatch: true)

        mutationsFiles =
            annotateMutation(byBamMutationChannel)
                .map { sampleId, file -> file }
                .collect()

        createMutationsTable(mutationsFiles, tumourMutationsChannel)

        offTargetErrorRates(createMutationsTable.out.filteredMutationsFile, layoutChannel)

        createOnTargetMutationsTable(
            createMutationsTable.out.filteredMutationsFile,
            tumourMutationsChannel,
            layoutChannel,
            offTargetErrorRates.out.noCosmicErrorRates)

        onTargetErrorRatesAndFilter(createOnTargetMutationsTable.out.onTargetMutationsFile,
                                    tumourMutationsChannel)
    emit:
        onTargetMutationsFile = onTargetErrorRatesAndFilter.out.onTargetMutationsFile
        onTargetLocusErrorRatesFile = onTargetErrorRatesAndFilter.out.locusErrorRates
        backgroundErrorRatesFile = createOnTargetMutationsTable.out.backgroundErrorRates
        offTargetErrorRatesWithCosmic = offTargetErrorRates.out.cosmicErrorRates
        offTargetErrorRatesNoCosmic = offTargetErrorRates.out.noCosmicErrorRates
}
