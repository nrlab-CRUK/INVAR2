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
        make_option(c("--alt-alleles-threashold"), type="integer", metavar="int",
                    dest="ALT_ALLELES_THRESHOLD", help="Blacklist loci with >= N separate alternate alleles",
                    default=3),
        make_option(c("--minor-alt-allele-threshold"), type="integer", metavar="int",
                    dest="MINOR_ALT_ALLELES_THRESHOLD", help="Blacklist multiallelic loci with a mutant read count of >= N in the minor mutant allele",
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
        ALT_ALLELES_THRESHOLD = 3,
        MINOR_ALT_ALLELES_THRESHOLD = 2
    )
}


##
# Loading functions.
#

# Load the patient specific tumour mutations file

loadTumourMutationsTable <- function(tumourMutationsFile)
{
    read_csv(tumourMutationsFile, col_types = 'ciccd', show_col_types = FALSE) %>%
    select(-contains('uniq'), -any_of('mut')) %>%
    rename_with(str_to_upper) %>%
    rename(CHROM = CHR) %>%
    mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':'))
}

# Read the layout file and extract unique pool id and barcode pairs.
# For this script, that column is all that is needed.

loadLayoutFile <- function(layoutFile)
{
    suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE)) %>%
        filter(case_or_control == "case") %>%
        mutate(POOL_BARCODE = str_c(SLX_ID, str_replace(barcode, '-', '_'), sep = '_')) %>%
        distinct(POOL_BARCODE)
}

# Load the mutations table from file.

loadMutationsTable <- function(mutationsFile, tumourMutationsTable, tapasSetting, cosmicThreshold)
{
    stopifnot(is.numeric(cosmicThreshold))

    # Remove soft-masked repeats, identified by lowercase.
    # Add columns of derived values and combined identifiers.
    read_tsv(mutationsFile, col_types = 'ciccicdddddccildc') %>%
        filter(!(str_detect(REF, '[acgt]') | str_detect(ALT, '[acgt]'))) %>%
        addDerivedColumns() %>%
        mutate(AF = (ALT_F + ALT_R) / DP,
               COSMIC = COSMIC_MUTATIONS > cosmicThreshold,
               SNP = `1KG_AF` > 0,
               ON_TARGET = UNIQUE_POS %in% tumourMutationsTable$UNIQUE_POS)
}

addDerivedColumns <- function(mutationTable)
{
    mutationTable %>%
        mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':'),
               POOL_BARCODE = str_c(POOL, BARCODE, sep='_'))
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
                                        minor_ALT_ALLELES_THRESHOLD)
{
    stopifnot(is.numeric(n_alt_alleles_threshold))
    stopifnot(is.numeric(minor_ALT_ALLELES_THRESHOLD))

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
        summarise(MUT_SUM = sum(ALT_F + ALT_R), .groups = "drop") %>%
        group_by(UNIQUE_POS) %>%
        filter(n() > 1) %>%
        summarise(N_ALT_ALLELES = n(),
                  MIN = min(MUT_SUM),
                  MAX = max(MUT_SUM),
                  .groups = "drop") %>%
        filter(N_ALT_ALLELES == n_alt_alleles_threshold &
               MIN >= minor_ALT_ALLELES_THRESHOLD)
}

##
# Saving functions
#

# Remove columns from the mutation table that can be derived from
# other columns, typically before saving.
removeDerivedColums <- function(mutationTable)
{
    mutationTable %>%
        select(-any_of(c('UNIQUE_POS', 'POOL_BARCODE')))
}

# Writes the TSV file but, before saving, converts any logical columns to
# simply the characters 'T' or 'F'. Saves having full "TRUE" and "FALSE" values,
# which are excessive as reading the table back correctly interprets 'T' and 'F'.
exportTSV <- function(t, file)
{
    toChar <- function(x)
    {
        ifelse(x, 'T', 'F')
    }

    t %>%
        mutate_if(is.logical, toChar) %>%
        write_tsv(file)

    t
}

# Save the given table as an RDS file and a TSV
saveRDSandTSV <- function(t, file)
{
    saveRDS(t, file)

    tsv <- str_replace(file, "\\.rds$", ".tsv")

    if (tsv != file)
    {
        exportTSV(t, tsv)
    }

    t
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    tumourMutationTable <- loadTumourMutationsTable(scriptArgs$TUMOUR_MUTATIONS_FILE)

    layoutTable <- loadLayoutFile(scriptArgs$LAYOUT_FILE)

    mutationTable.all <- loadMutationsTable(scriptArgs$MUTATIONS_FILE,
                                            tumourMutationTable,
                                            tapasSetting = scriptArgs$TAPAS_SETTING,
                                            cosmicThreshold = scriptArgs$COSMIC_THRESHOLD)

    # Filter the mutations table to remove rows with MSQB above the threshold
    # (was blacklist.MQSB in TAPAS_functions.R) and the base filters with
    # thresholds (IVAR3.R).

    mutationTable.filtered <- mutationTable.all %>%
        filter(MQSB > scriptArgs$MQSB_THRESHOLD &
               DP < scriptArgs$MAX_DP &
               !SNP &
               REF_R + REF_F >= scriptArgs$MIN_REF_DP &
               (ON_TARGET | !COSMIC))

    # Filter out positions that are multiallelic.

    multiallelicBlacklist <- mutationTable.filtered %>%
        createMultiallelicBlacklist(n_alt_alleles_threshold = scriptArgs$ALT_ALLELES_THRESHOLD,
                                    minor_ALT_ALLELES_THRESHOLD = scriptArgs$MINOR_ALT_ALLELES_THRESHOLD)

    mutationTable.biallelic <- mutationTable.filtered %>%
        filter(!UNIQUE_POS %in% multiallelicBlacklist$UNIQUE_POS)

    #mutationTable.all %>%
    #    removeDerivedColums() %>%
    #    saveRDSandTSV('mutation_table.all.rds')

    mutationTable.biallelic %>%
        removeDerivedColums() %>%
        saveRDSandTSV("mutation_table.filtered.rds")

    if (nrow(multiallelicBlacklist) > 0)
    {
        exportTSV(multiallelicBlacklist, 'multiallelic_blacklist.tsv')
    }
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789') {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
