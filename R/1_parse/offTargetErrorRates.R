suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))


##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"

    options_list <- list(
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by createMutationsTable.R",
                    default=defaultMarker),
        make_option(c("--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."),
        make_option(c("--control-proportion"), type="double", metavar="num",
                    dest="CONTROL_PROPORTION", help="Blacklist loci that have signal in >30% of the nonptspec samples",
                    default=0.1),
        make_option(c("--max-background-allele-frequency"), type="double", metavar="num",
                    dest="MAX_BACKGROUND_ALLELE_FREQUENCY", help="Filter loci with a background allele frequency in controls greater than this value",
                    default=0.01))

    opts <- OptionParser(option_list=options_list, usage="%prog [options] <mutation file>") %>%
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


##
# From TAPAS_functions.R, originally "annotate_with_locus_error_rate"
#
# Originally just annotate, this method creates the lociErrorRateTable and returns.
#
# Locus error rate = overall background error rate per locus, aggregated across control samples
#
createLociErrorRateTable <- function(mutationTable,
                                     layoutTable,
                                     proportionOfControls,
                                     maxBackgroundMeanAF,
                                     isBloodSpot)
{
    assert_that(is.number(proportionOfControls), msg = "proportionOfControls must be a number.")
    assert_that(is.number(maxBackgroundMeanAF), msg = "maxBackgroundMeanAF must be a number.")
    assert_that(is.flag(isBloodSpot), msg = "isBloodSpot must be a logical.")

    if (isBloodSpot)
    {
        message("sWGS/blood spot mutationTable, do not set a maxBackgroundMeanAF value as it is not appropriate in the low unique depth setting")
        maxBackgroundMeanAF <- 1
    }

    layoutTable.cases <- layoutTable %>%
        filter(CASE_OR_CONTROL == 'case')

    errorRateTable <- mutationTable %>%
        filter(SAMPLE_ID %in% layoutTable.cases$SAMPLE_ID) %>%
        mutate(HAS_SIGNAL = ifelse(ALT_F + ALT_R > 0, SAMPLE_ID, NA)) %>%
        group_by(CHROM, POS, TRINUCLEOTIDE) %>%
        summarize(MUTATION_SUM = sum(ALT_F) + sum(ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF = MUTATION_SUM / DP_SUM,
                  N_SAMPLES = n_distinct(SAMPLE_ID),
                  N_SAMPLES_WITH_SIGNAL = n_distinct(HAS_SIGNAL, na.rm = TRUE),
                  .groups = 'drop') %>%
        select(CHROM, POS, TRINUCLEOTIDE, BACKGROUND_AF, MUTATION_SUM, DP_SUM,
               N_SAMPLES, N_SAMPLES_WITH_SIGNAL) %>%
        mutate(LOCUS_NOISE.PASS = (N_SAMPLES_WITH_SIGNAL / N_SAMPLES) < proportionOfControls &
                                  BACKGROUND_AF < maxBackgroundMeanAF)

    errorRateTable
}

# Common function. Takes in a mutation table.
# Adds MUTATION_SUM, DP_SUM, BACKGROUND_AF columns on rows grouped by
# SAMPLE_ID, REF, ALT, TRINUCLEOTIDE
groupAndSummarizeForErrorRate <- function(mutationTable)
{
    mutationTable %>%
        group_by(SAMPLE_ID, REF, ALT, TRINUCLEOTIDE) %>%
        summarize(MUTATION_SUM = sum(ALT_F) + sum(ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF = (sum(ALT_F) + sum(ALT_R)) / sum(DP),
                  .groups = 'drop')

}

# Filter for rows that are off target, optionally also those that have COSMIC values.

filterForOffTarget <- function(mutationTable, withCosmic)
{
    assert_that(is.flag(withCosmic), msg = "withCosmic must be a logical.")

    mutationTable.off_target <- mutationTable %>%
        filter(!ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)

    if (!withCosmic)
    {
        mutationTable.off_target <- mutationTable.off_target %>%
            filter(!COSMIC)
    }

    mutationTable.off_target
}


# Add locus error rate and both strand information
addLocusNoisePass <- function(mutationTable, errorRateTable)
{
    passed.loci <- errorRateTable %>%
        filter(LOCUS_NOISE.PASS) %>%
        mutate(UNIQUE_POS = str_c(CHROM, POS, sep = ':'))

    mutationTable %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% passed.loci$UNIQUE_POS,
               BOTH_STRANDS.PASS = (ALT_R > 0 & ALT_F > 0) | AF == 0)
}


##
# Part of main, where this code is run with and without cosmic.
#

doMain <- function(withCosmic, mutationTable, layoutTable, lociErrorRateTable)
{
    assert_that(is.flag(withCosmic), msg = "withCosmic must be a logical.")

    cosmicFilePart = ifelse(withCosmic, 'cosmic', 'no_cosmic')

    mutationTable.off_target <- mutationTable %>%
        filterForOffTarget(withCosmic) %>%
        addLocusNoisePass(lociErrorRateTable)

    oneRead <- mutationTable.off_target %>%
        groupAndSummarizeForErrorRate()

    ## Calculate error rate with different settings and save

    # locus noise pass

    locusNoisePass <- mutationTable.off_target %>%
        filter(LOCUS_NOISE.PASS) %>%
        groupAndSummarizeForErrorRate()

    # both strands

    bothStrandsPass <- mutationTable.off_target %>%
        filter(BOTH_STRANDS.PASS) %>%
        groupAndSummarizeForErrorRate()

    # locus noise AND both strands

    bothFilters <- mutationTable.off_target %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS) %>%
        groupAndSummarizeForErrorRate()

    # Save these as a combined RDS file in a list.

    allErrorRates = list(
        PREFILTER = oneRead,
        LOCUS_NOISE = locusNoisePass,
        BOTH_STRANDS = bothStrandsPass,
        LOCUS_NOISE.BOTH_STRANDS = bothFilters
    )

    saveRDS(allErrorRates, str_c('error_rates.off_target.', cosmicFilePart, '.rds'))

    allErrorRates
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    assert_that(file.exists(scriptArgs$LAYOUT_FILE), msg = str_c(scriptArgs$LAYOUT_FILE, " does not exist."))
    assert_that(file.exists(scriptArgs$MUTATIONS_TABLE_FILE), msg = str_c(scriptArgs$MUTATIONS_TABLE_FILE, " does not exist."))

    layoutTable <-
        loadLayoutTable(scriptArgs$LAYOUT_FILE) %>%
        select(CASE_OR_CONTROL, SAMPLE_ID)

    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    # Doesn't matter about the COSMIC setting for the loci error rate table (no difference),
    # but the other filtering is necessary.

    lociErrorRateTable <- mutationTable %>%
        filterForOffTarget(FALSE) %>%
        createLociErrorRateTable(layoutTable,
                                 proportionOfControls = scriptArgs$CONTROL_PROPORTION,
                                 maxBackgroundMeanAF = scriptArgs$MAX_BACKGROUND_ALLELE_FREQUENCY,
                                 isBloodSpot = scriptArgs$BLOODSPOT)

    saveRDS(lociErrorRateTable, 'locus_error_rates.off_target.rds')

    # To make it as similar as possible to the old pipeline output when saved.
    # Aids comparison.

    lociErrorRateTable %>%
        select(-LOCUS_NOISE.PASS) %>%
        arrange(CHROM, POS, TRINUCLEOTIDE) %>%
        exportTSV('locus_error_rates.off_target.tsv')

    # Calculate the error rates with filters and with or without COSMIC.

    mclapply(c(TRUE, FALSE), doMain, mutationTable, layoutTable, lociErrorRateTable)
}

# Launch it.

invisible(main(parseOptions()))
