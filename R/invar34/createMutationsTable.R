suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))


##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"

    options_list <- list(
        make_option(c("--tapas"), type="character", metavar="string",
                    dest="TAPAS_SETTING", help="The TAPAS setting",
                    default=defaultMarker),
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_FILE", help="The calculated mutations file",
                    default=defaultMarker),
        make_option(c("--tumour-mutations"), type="character", metavar="file",
                    dest="TUMOUR_MUTATIONS_FILE", help="The source patient mutations file",
                    default=defaultMarker),
        make_option(c("--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("--cosmic-threshold"), type="integer", metavar="int",
                    dest="COSMIC_THRESHOLD", help="Loci with >0 entries in COSMIC are considered as COSMIC mutations",
                    default=0),
        make_option(c("--mqsb-threshold"), type="double", metavar="num",
                    dest="MQSB_THRESHOLD", help="Excludes data points due to poor MQ and SB, but locus is retained",
                    default=0.01),
        make_option(c("--max-dp"), type="integer", metavar="int",
                    dest="MAX_DP", help="Omit data points with uncharacteristially high unique DP given the input mass used",
                    default=1500),
        make_option(c("--min-ref-dp"), type="integer", metavar="int",
                    dest="MIN_REF_DP", help="Min DP is set by mpileups as 5, here we require at least 5 ref reads at a locus, set to 0 for sWGS",
                    default=5),
        make_option(c("--alt-allele-threashold"), type="integer", metavar="int",
                    dest="ALT_ALLELE_THRESHOLD", help="Blacklist loci with >= N separate alternate alleles",
                    default=3),
        make_option(c("--minor-alt-allele-threshold"), type="integer", metavar="int",
                    dest="MINOR_ALT_ALLELE_THRESHOLD", help="Blacklist multiallelic loci with a mutant read count of >= N in the minor mutant allele",
                    default=2))

    opts <- OptionParser(option_list=options_list, usage="%prog [options]") %>%
        parse_args(positional_arguments = TRUE)

    missingRequired <- which(opts$options == defaultMarker)
    if (any(missingRequired))
    {
        missing <- names(opts$options)[missingRequired]
        stop(str_c("ERROR: One or more required arguments have not been provided: ", str_c(missing, collapse = ' ')))
    }

    scriptOptions <- list()
    for (v in names(opts$options))
    {
        scriptOptions[v] = opts$options[v]
    }

    scriptOptions
}

# Test options for my (Rich) local set up in RStudio.

richTestOptions <- function()
{
    list(
        MUTATIONS_FILE = 'EMMA/output_gz/PARADIGM.f0.9_s2.BQ_20.MQ_40.combined.final.ann.tsv',
        TAPAS_SETTING = 'f0.9_s2.BQ_20.MQ_40',
        TUMOUR_MUTATIONS_FILE = 'source_files/PARADIGM_mutation_list_full_cohort_hg19.csv',
        LAYOUT_FILE = 'source_files/combined.SLX_table_with_controls_031220.csv',
        COSMIC_THRESHOLD = 0,
        MQSB_THRESHOLD = 0.01,
        MAX_DP = 2000,
        MIN_REF_DP = 10,
        ALT_ALLELE_THRESHOLD = 3,
        MINOR_ALT_ALLELE_THRESHOLD = 2
    )
}


##
# Taken from TAPAS_functions.R
# Multiallelic blacklisting.
#

# Filter out rows where the MSQB value is greater than the threshold.

blacklist.MQSB <- function(mutationTable, individual_MQSB_threshold)
{
    mutationTable %>%
        filter(MQSB > individual_MQSB_threshold)
}

# Create a list of loci that have more than the given number of alleles.

createMultiallelicBlacklist <- function(mutationTable, n_alt_alleles_threshold, minor_alt_allele_threshold)
{
    # Filter for rows with positive AF and
    # determine the number of alt alleles per UNIQUE_POS
    multiallelic <- mutationTable %>%
        filter(AF > 0) %>%
        mutate(MUT_SUM = ALT_F + ALT_R) %>%
        select(UNIQUE_POS, ALT, MUT_SUM) %>%
        arrange(UNIQUE_POS)

    # Identify loci with >1 alt alleles (e.g. both A>C and A>T at a position)
    duplicatePositions <- multiallelic %>%
        filter(duplicated(UNIQUE_POS))

    multiallelic <- multiallelic %>%
        filter(UNIQUE_POS %in% duplicatePositions)

    # Create a list of loci that have more than one allele.

    multiallelicBlacklist <- tibble(UNIQUE_POS = character(), .rows = 0)

    if (nrow(multiallelic) > 0)
    {
        multiallelicBlacklist <- multiallelic %>%
            group_by(UNIQUE_POS) %>%
            mutate(N_ALT_ALLELES = n(), MIN = min(MUT_SUM), MAX = max(MUT_SUM)) %>%
            filter(N_ALT_ALLELES >= n_alt_alleles_threshold &
                   MIN >= minor_alt_allele_threshold &
                   MAX >= minor_alt_allele_threshold) %>%
            select(UNIQUE_POS) %>%
            distinct(UNIQUE_POS)
    }

    multiallelicBlacklist
}

# Tidy multiallelic sites. Filter out those that have more than the requested
# number of alleles.

blacklist.multiallelic <- function(mutationTable,
                                   multiallelicBlacklist,
                                   n_alt_alleles_threshold,
                                   minor_alt_allele_threshold)
{
    filteredMutationTable <- mutationTable

    if (nrow(multiallelicBlacklist) > 0)
    {
        filteredMutationTable <- mutationTable %>%
            filter(!UNIQUE_POS %in% multiallelicBlacklist)
    }

    filteredMutationTable
}


##
# Loading functions.
#


# Load the patient specific tumour mutations file

load.tumour.mutations.table <- function(tumourMutationsFile)
{
    read_csv(tumourMutationsFile, col_types = 'ciccd', show_col_types = FALSE) %>%
    select(-contains('uniq'), -mut) %>%
    rename_with(str_to_upper) %>%
    rename(CHROM = CHR) %>%
    mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':'))
}

# Read the layout file and extract unique pool id and barcode pairs.

load.layout.file <- function(layoutFile)
{
    suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE)) %>%
        filter(case_or_control == "case") %>%
        mutate(POOL_BARCODE = str_c(SLX_ID, str_replace(barcode, '-', '_'), sep = '_')) %>%
        select(POOL_BARCODE) %>%
        distinct(POOL_BARCODE)
}

addDerivedColumns <- function(mutationTable)
{
    mutationTable %>%
        mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':')) %>%
        mutate(POOL_BARCODE = str_c(POOL, BARCODE, sep='_')) %>%
        mutate(UNIQUE_ID = str_c(UNIQUE_POS, POOL_BARCODE, sep='_'))
}

# Load the mutations table from file.

load.mutations.table <- function(mutationsFile, tumourMutationsTable, tapasSetting, cosmic_threshold)
{
    message("Reading mutations table from ", mutationsFile)

    # Remove soft-masked repeats, identified by lowercase.
    # Add columns of derived values and combined identifiers.
    read_tsv(mutationsFile, col_types = 'ciccicdddddccildc') %>%
        filter(!(str_detect(REF, '[acgt]') | str_detect(ALT, '[acgt]'))) %>%
        addDerivedColumns() %>%
        mutate(AF = (ALT_F + ALT_R) / DP) %>%
        mutate(COSMIC = COSMIC_MUTATIONS > cosmic_threshold) %>%
        mutate(SNP = `1KG_AF` > 0) %>%
        mutate(ON_TARGET = UNIQUE_POS %in% tumourMutationsTable$UNIQUE_POS)
}

# Filter the mutations table to remove rows with MSQB above the threshold
# and positions that are multiallelic.

filter.mutations.table <- function(mutationTable,
                                   tumourMutationsTable,
                                   multiallelicBlacklist,
                                   individual_MQSB_threshold,
                                   max_DP,
                                   min_ref_DP,
                                   n_alt_alleles_threshold,
                                   minor_alt_allele_threshold)
{
    # Apply filters to blacklist loci
    mutationTable %>%
        filter(DP < max_DP & !SNP & REF_R + REF_F >= min_ref_DP & (ON_TARGET | !COSMIC)) %>%
        blacklist.MQSB(individual_MQSB_threshold) %>%
        blacklist.multiallelic(multiallelicBlacklist, n_alt_alleles_threshold, minor_alt_allele_threshold)
}

saveRDSandTSV <- function(t, file)
{
    cleaned <- t %>%
        select(-any_of(c('UNIQUE_POS', 'POOL_BARCODE', 'UNIQUE_ID')))

    saveRDS(cleaned, file)

    tsv <- str_replace(file, "\\.rds$", ".tsv")

    if (tsv != file)
    {
        write_tsv(cleaned, tsv)
    }

    t
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    tumourMutationTable <- load.tumour.mutations.table(scriptArgs$TUMOUR_MUTATIONS_FILE)

    layoutTable <- load.layout.file(scriptArgs$LAYOUT_FILE)

    mutationTable.all <- load.mutations.table(scriptArgs$MUTATIONS_FILE,
                                              tumourMutationTable,
                                              scriptArgs$TAPAS_SETTING,
                                              scriptArgs$COSMIC_THRESHOLD)

    saveRDSandTSV(mutationTable.all, 'mutation_table.all.rds')

    multiallelicBlacklist <- mutationTable.all %>%
        blacklist.MQSB(scriptArgs$MQSB_THRESHOLD) %>%
        createMultiallelicBlacklist(scriptArgs$ALT_ALLELES_THRESHOLD, scriptArgs$MINOR_ALT_ALLELE_THRESHOLD)

    mutationTable.all %>%
        filter.mutations.table(tumourMutationTable,
                               multiallelicBlacklist,
                               individual_MQSB_threshold = scriptArgs$MQSB_THRESHOLD,
                               max_DP = scriptArgs$MAX_DP,
                               min_ref_DP = scriptArgs$MIN_REF_DP,
                               n_alt_alleles_threshold = scriptArgs$ALT_ALLELES_THRESHOLD,
                               minor_alt_allele_threshold = scriptArgs$MINOR_ALT_ALLELE_THRESHOLD) %>%
        saveRDSandTSV("mutation_table.filtered.rds")

    if (nrow(multiallelicBlacklist) > 0)
    {
        write_tsv(multiallelicBlacklist, 'multiallelic_blacklist.tsv', col_names = FALSE)
    }
}

# Launch it.

invisible(main(parseOptions()))
#invisible(main(richTestOptions()))
