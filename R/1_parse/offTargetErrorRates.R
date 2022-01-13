suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/1_parse/common.R'))


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

# Test options for my (Rich) local set up in RStudio.

richTestOptions <- function()
{
    list(
        MUTATIONS_TABLE_FILE = 'mutations/mutation_table.filtered.rds',
        LAYOUT_FILE = 'source_files/combined.SLX_table_with_controls_031220.csv',
        CONTROL_PROPORTION = 0.1,
        MAX_BACKGROUND_ALLELE_FREQUENCY = 0.01,
        BLOODSPOT = FALSE
    )
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
                                     proportion_of_controls,
                                     max_background_mean_AF,
                                     is.blood_spot)
{
    assert_that(is.number(proportion_of_controls), msg = "proportion_of_controls must be a number.")
    assert_that(is.number(max_background_mean_AF), msg = "max_background_mean_AF must be a number.")
    assert_that(is.flag(is.blood_spot), msg = "is.blood_spot must be a logical.")

    if (is.blood_spot)
    {
        message("sWGS/blood spot mutationTable, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
        max_background_mean_AF <- 1
    }

    layoutTable.cases <- layoutTable %>%
        filter(CASE_OR_CONTROL == 'case')

    errorRateTable <- mutationTable %>%
        filter(POOL_BARCODE %in% layoutTable.cases$POOL_BARCODE) %>%
        mutate(HAS_SIGNAL = ifelse(ALT_F + ALT_R > 0, POOL_BARCODE, NA)) %>%
        group_by(UNIQUE_POS, TRINUCLEOTIDE) %>%
        summarize(MUTATION_SUM = sum(ALT_F) + sum(ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF = MUTATION_SUM / DP_SUM,
                  N_SAMPLES = n_distinct(POOL_BARCODE),
                  N_SAMPLES_WITH_SIGNAL = n_distinct(HAS_SIGNAL, na.rm = TRUE),
                  .groups = 'drop') %>%
        mutate(LOCUS_NOISE.PASS = (N_SAMPLES_WITH_SIGNAL / N_SAMPLES) < proportion_of_controls &
                                  BACKGROUND_AF < max_background_mean_AF)

    errorRateTable
}

# Common function. Takes in a mutation table.
# Adds MUTATION_SUM, DP_SUM, BACKGROUND_AF columns on rows grouped by
# POOL, BARCODE, REF, ALT, TRINUCLEOTIDE
groupAndSummarizeForErrorRate <- function(mutationTable)
{
    mutationTable %>%
        group_by(POOL, BARCODE, REF, ALT, TRINUCLEOTIDE) %>%
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
        filter(LOCUS_NOISE.PASS)

    mutationTable %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% passed.loci$UNIQUE_POS,
               BOTH_STRANDS = (ALT_R > 0 & ALT_F > 0) | AF == 0)
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

    bothStrands <- mutationTable.off_target %>%
        filter(BOTH_STRANDS) %>%
        groupAndSummarizeForErrorRate()

    # locus noise AND both strands

    bothFilters <- mutationTable.off_target %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS) %>%
        groupAndSummarizeForErrorRate()

    # Save these as a combined RDS file in a list.

    allErrorRates = list(
        ONE_READ = oneRead,
        LOCUS_NOISE = locusNoisePass,
        BOTH_READS = bothStrands,
        LOCUS_NOISE.BOTH_READS = bothFilters
    )

    saveRDS(allErrorRates, str_c('mutation_table.error_rates.', cosmicFilePart, '.rds'))

    # Save tables as TSV for reference

    if (TRUE)
    {
        oneRead %>%
            removeMutationTableDerivedColumns() %>%
            exportTSV(str_c('mutation_table.off_target.', cosmicFilePart, '.oneread.tsv'))

        locusNoisePass %>%
            removeMutationTableDerivedColumns() %>%
            exportTSV(str_c('mutation_table.off_target.', cosmicFilePart, '.locusnoise.tsv'))

        bothStrands %>%
            removeMutationTableDerivedColumns() %>%
            exportTSV(str_c('mutation_table.off_target.', cosmicFilePart, '.bothreads.tsv'))

        bothFilters %>%
            removeMutationTableDerivedColumns() %>%
            exportTSV(str_c('mutation_table.off_target.', cosmicFilePart, '.locusnoise_bothreads.tsv'))
    }

    allErrorRates
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    layoutTable <-
        loadLayoutTable(scriptArgs$LAYOUT_FILE) %>%
        filter(CASE_OR_CONTROL == "case")

    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    # Doesn't matter about the COSMIC setting for the loci error rate table (no difference),
    # but the other filtering is necessary.

    lociErrorRateTable <- mutationTable %>%
        filterForOffTarget(FALSE) %>%
        createLociErrorRateTable(layoutTable,
                                 proportion_of_controls = scriptArgs$CONTROL_PROPORTION,
                                 max_background_mean_AF = scriptArgs$MAX_BACKGROUND_ALLELE_FREQUENCY,
                                 is.blood_spot = scriptArgs$BLOODSPOT)

    saveRDS(lociErrorRateTable, 'locus_error_rates.off_target.rds')

    # To make it as similar as possible to the old pipeline output when saved.
    # Aids comparison.

    lociErrorRateTable %>%
        select(UNIQUE_POS, TRINUCLEOTIDE, BACKGROUND_AF, MUTATION_SUM, DP_SUM, N_SAMPLES, N_SAMPLES_WITH_SIGNAL) %>%
        exportTSV('locus_error_rates.off_target.tsv')

    # Calculate the error rates with filters and with or without COSMIC.

    mclapply(c(TRUE, FALSE), doMain, mutationTable, layoutTable, lociErrorRateTable)
    #doMain(TRUE, mutationTable, layoutTable, lociErrorRateTable)
    #doMain(FALSE, mutationTable, layoutTable, lociErrorRateTable)
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789') {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
