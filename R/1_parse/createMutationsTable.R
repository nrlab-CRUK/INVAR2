suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))


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
                    default=0L),
        make_option(c("--mqsb-threshold"), type="double", metavar="num",
                    dest="MQSB_THRESHOLD", help="Excludes data points due to poor MQ and SB, but locus is retained",
                    default=0.01),
        make_option(c("--max-depth"), type="integer", metavar="int",
                    dest="MAX_DEPTH", help="Omit data points with uncharacteristially high unique depth given the input mass used",
                    default=1500L),
        make_option(c("--min-ref-depth"), type="integer", metavar="int",
                    dest="MIN_REF_DEPTH", help="Min depth is set by mpileups as 5, here we require at least 5 ref reads at a locus, set to 0 for sWGS",
                    default=5L),
        make_option(c("--alt-alleles-threshold"), type="integer", metavar="int",
                    dest="ALT_ALLELES_THRESHOLD", help="Blacklist loci with >= N separate alternate alleles",
                    default=3L),
        make_option(c("--minor-alt-allele-threshold"), type="integer", metavar="int",
                    dest="MINOR_ALT_ALLELES_THRESHOLD", help="Blacklist multiallelic loci with a mutant read count of >= N in the minor mutant allele",
                    default=2L))

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
        COSMIC_THRESHOLD = 0L,
        MQSB_THRESHOLD = 0.01,
        MAX_DEPTH = 2000L,
        MIN_REF_DEPTH = 10L,
        ALT_ALLELES_THRESHOLD = 3L,
        MINOR_ALT_ALLELES_THRESHOLD = 2L
    )
}


##
# Loading functions.
#

# Load the mutations table from file.

loadMutationsTable <- function(mutationsFile, tumourMutationsTable, tapasSetting, cosmicThreshold)
{
    assert_that(is.number(cosmicThreshold), msg = "cosmicThreshold must be a number.")

    # Remove soft-masked repeats, identified by lowercase.
    # Add columns of derived values and combined identifiers.
    read_tsv(mutationsFile, col_types = 'cicciciiiidccildc') %>%
        filter(!(str_detect(REF, '[acgt]') | str_detect(ALT, '[acgt]'))) %>%
        addMutationTableDerivedColumns() %>%
        mutate(AF = (ALT_F + ALT_R) / DP,
               COSMIC = COSMIC_MUTATIONS > cosmicThreshold,
               SNP = `1KG_AF` > 0,
               ON_TARGET = UNIQUE_POS %in% tumourMutationsTable$UNIQUE_POS)
}

##
# Filtering functions.
#

##
# Taken from TAPAS_functions.R
# Multiallelic blacklisting.
#
# Create a list of loci that have more than the given number of alleles.

createMultiallelicBlacklist <- function(mutationTable,
                                        n_alt_alleles_threshold,
                                        minor_alt_alleles_threshold)
{
    assert_that(is.number(n_alt_alleles_threshold), msg = "n_alt_alleles_threshold must be a number")
    assert_that(is.number(minor_alt_alleles_threshold), msg = "minor_alt_alleles_threshold must be a number")

    # Filter for rows with positive AF and MQSB above threshold.
    # Then determine the number of alt alleles per UNIQUE_POS
    # Identify loci with >1 alt alleles (e.g. both A>C and A>T at a position)

    # In the original code, the final filter is:
    # N_ALT_ALLELES == n_alt_alleles_threshold
    # I'm sure this should be >=, but to keep things the same for now,
    # it's == here too.

    mutationTable %>%
        filter(AF > 0) %>%
        group_by(UNIQUE_POS, ALT) %>%
        summarise(MUTATION_SUM = sum(ALT_F + ALT_R), .groups = "drop") %>%
        group_by(UNIQUE_POS) %>%
        filter(n() > 1) %>%
        summarise(N_ALT_ALLELES = n(),
                  MIN = min(MUTATION_SUM),
                  MAX = max(MUTATION_SUM),
                  .groups = "drop") %>%
        filter(N_ALT_ALLELES == n_alt_alleles_threshold &
               MIN >= minor_alt_alleles_threshold)
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    print(scriptArgs)

    tumourMutationTable <-
        loadTumourMutationsTable(scriptArgs$TUMOUR_MUTATIONS_FILE)

    # Read the layout file and extract unique pool id and barcode pairs.
    # For this script, that column is all that is needed.

    layoutTable <-
        loadLayoutTable(scriptArgs$LAYOUT_FILE) %>%
        distinct(POOL_BARCODE)

    mutationTable.all <- loadMutationsTable(scriptArgs$MUTATIONS_FILE,
                                            tumourMutationTable,
                                            tapasSetting = scriptArgs$TAPAS_SETTING,
                                            cosmicThreshold = scriptArgs$COSMIC_THRESHOLD)

    message("Number of mutations at start = ", nrow(mutationTable.all))

    # Filter the mutations table to remove rows with MSQB above the threshold
    # (was blacklist.MQSB in TAPAS_functions.R) and the base filters with
    # thresholds (INVAR3.R).

    mutationTable.filtered <- mutationTable.all %>%
        filter(!is.na(MQSB) &
               MQSB > scriptArgs$MQSB_THRESHOLD &
               DP < scriptArgs$MAX_DEPTH &
               !SNP &
               REF_R + REF_F >= scriptArgs$MIN_REF_DEPTH &
               (ON_TARGET | !COSMIC))

    message("Number of mutations after MQSB filter = ", nrow(mutationTable.filtered))

    # Filter out positions that are multiallelic.

    multiallelicBlacklist <- mutationTable.filtered %>%
        createMultiallelicBlacklist(n_alt_alleles_threshold = scriptArgs$ALT_ALLELES_THRESHOLD,
                                    minor_alt_alleles_threshold = scriptArgs$MINOR_ALT_ALLELES_THRESHOLD)

    message("Number of multiallelic blacklist positions = ", nrow(multiallelicBlacklist))

    mutationTable.biallelic <- mutationTable.filtered %>%
        filter(!UNIQUE_POS %in% multiallelicBlacklist$UNIQUE_POS)

    message("Number of mutations after multiallelic filter = ", nrow(mutationTable.biallelic))

    #mutationTable.all %>%
    #    removeMutationTableDerivedColumns() %>%
    #    arrange(POOL, BARCODE, CHROM, POS, REF, ALT, TRINUCLEOTIDE) %>%
    #    saveRDSandTSV('mutation_table.all.rds')

    mutationTable.biallelic %>%
        removeMutationTableDerivedColumns() %>%
        arrange(POOL, BARCODE, CHROM, POS, REF, ALT, TRINUCLEOTIDE) %>%
        saveRDSandTSV("mutation_table.filtered.rds")

    if (nrow(multiallelicBlacklist) > 0)
    {
        exportTSV(multiallelicBlacklist, 'multiallelic_blacklist.tsv')
    }
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine in RStudio
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
