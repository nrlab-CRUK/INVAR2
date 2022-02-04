process slopPatientInfo
{
    executor 'local'
    memory '256m'
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
    tag "${pool} ${barcode}"

    cpus { Math.min(params.MAX_CORES, 2) }

    input:
        each path(sloppedBedFile)
        each path(fastaReference)
        tuple val(pool), val(barcode), path(bamFile), path(bamIndex)

    output:
        tuple val(pool), val(barcode), path(vcfFile), emit: "vcfFile"

    shell:
        vcfFile = "${pool}.${barcode}.mpileup.vcf"

        dedupFlags = params.REMOVE_DUPLICATES ? "-R --ff UNMAP" : ""

        template "1_parse/mpileup.sh"
}


process biallelic
{
    tag "${pool} ${barcode}"

    memory '256m'
    time   '10m'

    input:
        tuple val(pool), val(barcode), path(vcfFile)

    output:
        tuple val(pool), val(barcode), path(mutationsFile), emit: "mutationsFile"

    shell:
        mutationsFile = "${pool}.${barcode}.tsv"

        template "1_parse/biallelic.sh"
}


process tabixSnp
{
    tag "${pool} ${barcode}"

    memory '32m'

    input:
        each path(tabixDatabase)
        each path(tabixDatabaseIndex)
        tuple val(pool), val(barcode), path(mutationFile)

    output:
        tuple val(pool), val(barcode), path(tabixFile)

    shell:
        tabixFile = "${pool}.${barcode}.snp.vcf"

        template "1_parse/tabix.sh"
}

process tabixCosmic
{
    tag "${pool} ${barcode}"

    memory '32m'

    input:
        each path(tabixDatabase)
        each path(tabixDatabaseIndex)
        tuple val(pool), val(barcode), path(mutationFile)

    output:
        tuple val(pool), val(barcode), path(tabixFile)

    shell:
        tabixFile = "${pool}.${barcode}.cosmic.vcf"

        template "1_parse/tabix.sh"
}

process trinucleotide
{
    tag "${pool} ${barcode}"

    memory '256m'

    input:
        each path(fastaReference)
        tuple val(pool), val(barcode), path(mutationFile)

    output:
        tuple val(pool), val(barcode), path(trinucleotideFile)

    shell:
        trinucleotideFile = "${pool}.${barcode}.trinucleotide.fa"

        template "1_parse/samtools_faidx.sh"
}

process annotateMutation
{
    tag "${pool} ${barcode}"

    input:
        tuple val(pool), val(barcode), path(mutationFile), path(snp), path(cosmic), path(trinucleotide)

    output:
        tuple val(pool), val(barcode), path(annotatedFile), emit: "mutationsFile"

    shell:
        annotatedFile = "${pool}.${barcode}.mutations.tsv"

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

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile

    output:
        path 'mutation_table.filtered.rds', emit: "filteredMutationsFile"

    shell:
        """
        Rscript --vanilla "!{params.projectHome}/R/1_parse/createMutationsTable.R" \
            --mutations="!{mutationsFile}" \
            --tumour-mutations="!{tumourMutationsFile}" \
            --layout="!{layoutFile}" \
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
    cpus   { Math.min(params.MAX_CORES, 2) }

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
            --control-proportion=!{params.proportion_of_controls} \
            --max-background-allele-frequency=!{params.max_background_mean_allele_frequency} \
            !{params.is_bloodspot ? "--bloodspot" : ""}
        """
}

process createOnTargetMutationsTable
{
    memory '4g'

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
    memory '4g'

    input:
        path mutationsFile
        path tumourMutationsFile
        path layoutFile

    output:
        path 'locus_error_rates.on_target.rds', emit: 'locusErrorRates'
        path 'locus_error_rates.on_target.pdf', optional: true, emit: 'locusErrorRatesPlot'
        path 'mutation_table.on_target.rds', emit: "onTargetMutationsFile"

    shell:
        tapasSetting = "${params.ERROR_SUPPRESSION_NAME}_BQ_${params.BASEQ}.MQ_${params.MAPQ}"

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
        genome_channel = channel.fromPath(params.HG19_GENOME, checkIfExists: true)
        fasta_channel = channel.fromPath(params.FASTAREF, checkIfExists: true)
        snp_channel = channel.fromPath(params.K1G_DB, checkIfExists: true)
        snp_index_channel = channel.fromPath("${params.K1G_DB}.tbi")
        cosmic_channel = channel.fromPath(params.COSMIC_DB, checkIfExists: true)
        cosmic_index_channel = channel.fromPath("${params.COSMIC_DB}.tbi")

        slopPatientInfo(tumourMutationsChannel, genome_channel)

        mpileup(slopPatientInfo.out, fasta_channel, bamChannel)
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
        onTargetErrorRatesFile = onTargetErrorRatesAndFilter.out.locusErrorRates
        backgroundErrorRatesFile = createOnTargetMutationsTable.out.backgroundErrorRates
}
