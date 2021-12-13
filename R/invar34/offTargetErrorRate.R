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
        make_option(c("-l", "--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("-b", "--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."),
        make_option(c("-c", "--cosmic"), action="store_true", default=FALSE,
                    dest="COSMIC", help="Keep Cosmic scores."))

    opts <- OptionParser(option_list=options_list, usage="%prog [options] <mutation file>") %>%
        parse_args(positional_arguments = TRUE)

    missingRequired <- which(opts$options == defaultMarker)
    if (any(missingRequired))
    {
        missing <- names(opts$options)[missingRequired]
        stop(str_c("ERROR: One or more required arguments have not been provided: ", str_c(missing, collapse = ' ')))
    }

    if (length(opts$args) != 1)
    {
        stop("The mutation table file (RDS) has not been provided.")
    }

    scriptOptions <- list(MUTATIONS_TABLE_FILE = opts$args[1])
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
        MUTATIONS_TABLE_FILE = 'invar34/mutation_table.filtered.rds',
        LAYOUT_FILE = 'source_files/combined.SLX_table_with_controls_031220.csv',
        BLOODSPOT = FALSE,
        COSMIC = TRUE
    )
}


##
# Loading functions.
#

# Read the layout file and extract unique pool id and barcode pairs.

load.layout.file <- function(layoutFile)
{
    suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE)) %>%
    filter(case_or_control == "case") %>%
    mutate(POOL_BARCODE = str_c(SLX_ID, str_replace(barcode, '-', '_'), sep = '_')) %>%
    select(POOL_BARCODE) %>%
    distinct(POOL_BARCODE)
}

##
# From TAPAS_functions.R, originally "annotate_with_locus_error_rate"
#
# Originally just annotate, this method creates the errorRateTable and returns.
#
# Locus error rate = overall background error rate per locus, aggregated across control samples
createErrorRateTable <- function(mutationTable,
                                 layoutTable,
                                 on_target,
                                 proportion_of_controls = 0.1,
                                 max_background_mean_AF = 0.01,
                                 is.blood_spot = FALSE)
{
    if (is.blood_spot)
    {
        message("sWGS/blood spot mutationTable, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
        max_background_mean_AF <- 1
    }

    errorRateTable <- mutationTable %>%
        filter(POOL_BARCODE %in% layoutTable$POOL_BARCODE) %>%
        mutate(HAS_SIGNAL = ifelse(ALT_F + ALT_R > 0, POOL_BARCODE, NA)) %>%
        group_by(UNIQUE_POS, TRINUCLEOTIDE) %>%
        summarize(MUT_SUM = sum(ALT_F) + sum(ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF = MUT_SUM / DP_SUM,
                  N_SAMPLES = n_distinct(POOL_BARCODE),
                  N_SAMPLES_WITH_SIGNAL = n_distinct(HAS_SIGNAL, na.rm = TRUE),
                  .groups = 'drop')

    errorRateTable %>%
        mutate(PROPORTION = N_SAMPLES_WITH_SIGNAL / N_SAMPLES,
               LOCUS_NOISE.PASS = PROPORTION < proportion_of_controls & BACKGROUND_AF < max_background_mean_AF,
               HAS_AF = BACKGROUND_AF > 0,
               LOCUS_NOISE.FAIL = !LOCUS_NOISE.PASS)
}

# Common function. Takes in a mutation table.
# Adds MUT_SUM, DP_SUM, BACKGROUND_AF columns on rows grouped by
# POOL, BARCODE, REF, ALT, TRINUCLEOTIDE
groupAndSummarizeForErrorRate <- function(mt)
{
    mt %>%
        group_by(POOL, BARCODE, REF, ALT, TRINUCLEOTIDE) %>%
        summarize(MUT_SUM = sum(ALT_F) + sum(ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF= (sum(ALT_F) + sum(ALT_R)) / sum(DP),
                  .groups = 'drop')

}

# Filter for rows that are off target, optionally also those that have COSMIC values.

filter.off_target <- function(mutationTable, use_cosmic)
{
    mutationTable.off_target <- mutationTable %>%
        filter(!ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)

    if (!use_cosmic)
    {
        mutationTable.off_target <- mutationTable.off_target %>%
            filter(!COSMIC)
    }

    mutationTable.off_target
}


# Add locus error rate and both strand information
add.locus_noise_pass <- function(mutationTable, errorRateTable)
{
    passed.loci <- errorRateTable %>%
        filter(LOCUS_NOISE.PASS) %>%
        select(UNIQUE_POS)

    mutationTable %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% passed.loci$UNIQUE_POS,
               BOTH_STRANDS = (ALT_R > 0 & ALT_F > 0) | AF == 0)
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    cosmicFilePart = ifelse(scriptArgs$COSMIC, 'cosmic', 'no_cosmic')

    layoutTable <- load.layout.file(scriptArgs$LAYOUT_FILE)

    mutationTable <- readRDS(scriptArgs$MUTATIONS_TABLE_FILE)

    errorRateTable <- createErrorRateTable(mutationTable, layoutTable, FALSE, scriptArgs$BLOODSPOT)

    mutationTable.off_target <- mutationTable %>%
        filter.off_target(scriptArgs$COSMIC) %>%
        add.locus_noise_pass(errorRateTable)

    mutationTable.off_target %>%
        groupAndSummarize() %>%
        writeRDS(str_c('locus_error_rates.off_target.', cosmicFilePart, '.all.rds'))

    ## Calculate error rate with different settings and save

    # locus noise pass

    mutations.off_target %>%
        filter(LOCUS_NOISE.PASS) %>%
        groupAndSummarize() %>%
        writeRDS(str_c('mutation_table.off_target.', cosmicFilePart, '.locusnoisepass.rds'))

    # both strands

    mutations.off_target %>%
        filter(BOTH_STRANDS) %>%
        groupAndSummarize() %>%
        writeRDS(str_c('mutation_table.off_target.', cosmicFilePart, '.bothstrands.rds'))

    # locus noise AND both strands

    background_error.locus_noise.both_strands <- mutations.off_target %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS) %>%
        groupAndSummarize() %>%
        writeRDS(str_c('mutation_table.off_target.', cosmicFilePart, '.locusnoisepass_bothstrands.rds'))
}

# Launch it.

invisible(main(parseOptions()))
