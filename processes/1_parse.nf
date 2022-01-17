// tapasSetting is a closure that can be used inside more than one process that
// gives the full string for the TAPAS setting.

tapasSetting =
{
    return "${params.ERROR_SUPPRESSION_NAME}_BQ_${params.BASEQ}.MQ_${params.MAPQ}"
}

// Processes

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
        python3 "!{params.projectHome}/python/1_parse/addTabixAndTrinucleotides.py" \
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
        combinedFile = "mutation_table.tsv"

        template "1_parse/catCSV.sh"
}

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
        path 'mutation_table.filtered.tsv', emit: "filteredMutationsTSV"

    shell:

        """
        Rscript --vanilla "!{params.projectHome}/R/1_parse/createMutationsTable.R" \
            --mutations="!{mutationsFile}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
            --tapas="!{tapasSetting}" \
            --cosmic-threshold=!{params.cosmic_threshold} \
            --mqsb-threshold=!{params.individual_MQSB_threshold} \
            --max-depth=!{params.max_depth} \
            --min-ref-depth=!{params.min_ref_depth} \
            --alt-alleles-threshold=!{params.n_alt_alleles_threshold} \
            --minor-alt-allele-threshold=!{params.minor_alt_allele_threshold}
        """
}

process offTargetErrorRates
{
    memory '4g'
    cpus 2
    time '1h'

    publishDir 'off_target', mode: 'link'

    input:
        path mutationsFile
        path layoutFile

    output:
        path 'locus_error_rates.off_target.rds', emit: 'locusErrorRates'
        path 'locus_error_rates.off_target.tsv', emit: 'locusErrorRatesTSV'
        path 'mutation_table.error_rates.cosmic.rds', emit: "cosmicErrorRates"
        path 'mutation_table.error_rates.no_cosmic.rds', emit: "noCosmicErrorRates"
        path 'mutation_table.off_target.*.tsv', emit: "errorRatesTSV"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/1_parse/offTargetErrorRates.R" \
            --mutations="!{mutationsFile}" \
            --layout="!{layoutFile}" \
            --control-proportion=!{params.proportion_of_controls} \
            --max-background-allele-frequency=!{params.max_background_mean_allele_frequency} \
            !{params.is_bloodspot ? "--bloodspot" : ""}
        """
}

process createOnTargetMutationsTable
{
    memory '4g'
    cpus 1
    time '1h'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile
        path errorRatesFile

    output:
        path 'mutation_table.on_target.all.rds', emit: "onTargetMutationsFile"
        path '*.tsv', optional: true

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
    memory '4g'
    cpus 1
    time '1h'

    publishDir 'on_target', mode: 'link'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile

    output:
        path 'locus_error_rates.on_target.rds', emit: 'locusErrorRates'
        path 'mutation_table.on_target.rds', emit: "onTargetMutationsFile"
        path '*.tsv', optional: true
        path '*.pdf', optional: true

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/1_parse/onTargetErrorRatesAndFilter.R" \
            --mutations="!{mutationsFile}" \
            --layout="!{layoutFile}" \
            --study="!{params.STUDY_ID}" \
            --tapas="!{tapasSetting}" \
            --cosmic-threshold=!{params.cosmic_threshold} \
            --control-proportion=!{params.proportion_of_controls} \
            --max-background-allele-frequency=!{params.max_background_mean_allele_frequency} \
            !{params.is_bloodspot ? "--bloodspot" : ""} \
            --allele-frequency-threshold=!{params.allele_frequency_threshold}
        """
}


workflow parse
{
    take:
        bamChannel
        tumourMutationsChannel
        layoutChannel

    main:
        simpleBamChannel = bamChannel.map { pool, barcode, bam, index -> bam }

        genome_channel = channel.fromPath(params.HG19_GENOME, checkIfExists: true)
        fasta_channel = channel.fromPath(params.FASTAREF, checkIfExists: true)
        snp_channel = channel.fromPath(params.K1G_DB, checkIfExists: true)
        snp_index_channel = channel.fromPath("${params.K1G_DB}.tbi")
        cosmic_channel = channel.fromPath(params.COSMIC_DB, checkIfExists: true)
        cosmic_index_channel = channel.fromPath("${params.COSMIC_DB}.tbi")

        slopPatientInfo(tumourMutationsChannel, genome_channel)

        mpileup(slopPatientInfo.out, fasta_channel, simpleBamChannel)
        biallelic(mpileup.out)

        mutationChannel = biallelic.out.filter { pool, bc, f -> f.countLines() > 1 }

        tabixSnp(snp_channel, snp_index_channel, mutationChannel)
        tabixCosmic(cosmic_channel, cosmic_index_channel, mutationChannel)
        trinucleotide(fasta_channel, mutationChannel)

        byBamMutationChannel = mutationChannel
            .join(tabixSnp.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)
            .join(tabixCosmic.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)
            .join(trinucleotide.out, by: 0..1, failOnDuplicate: true, failOnMismatch: true)

        allMutationsList =
            annotateMutation(byBamMutationChannel)
                .map { pool, barcode, file -> file }
                .toSortedList( { f1, f2 -> f1.name <=> f2.name } )

        combineCSV(allMutationsList)

        createMutationsTable(combineCSV.out, tumourMutationsChannel, layoutChannel)

        offTargetErrorRates(createMutationsTable.out.filteredMutationsFile, layoutChannel)

        createOnTargetMutationsTable(
            createMutationsTable.out.filteredMutationsFile,
            tumourMutationsChannel,
            layoutChannel,
            offTargetErrorRates.out.noCosmicErrorRates)

        onTargetErrorRatesAndFilter(createOnTargetMutationsTable.out.onTargetMutationsFile,
                                    tumourMutationsChannel,
                                    layoutChannel)
    emit:
        mutationsFile = combineCSV.out
        onTargetMutationsFile = onTargetErrorRatesAndFilter.out.onTargetMutationsFile
}
