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
        make_option(c("-t", "--tapas"), type="character", metavar="string",
                    dest="TAPAS_SETTING", help="The TAPAS setting",
                    default=defaultMarker),
        make_option(c("-m", "--tumour-mutations"), type="character", metavar="file",
                    dest="TUMOUR_MUTATIONS_FILE", help="The source patient mutations file",
                    default=defaultMarker),
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
    #mutationsFile <- 'mutations/SLX-19722_SXTLI005.mutations.tsv'
    list(
        MUTATIONS_FILE = '/mnt/scratchb/bioinformatics/bowers01/invar_emma/EMMA/output_gz/PARADIGM.f0.9_s2.BQ_20.MQ_40.combined.final.ann.tsv',
        TAPAS_SETTING = 'f0.9_s2.BQ_20.MQ_40',
        TUMOUR_MUTATIONS_FILE = 'source_files/PARADIGM_mutation_list_full_cohort_hg19.csv',
        LAYOUT_FILE = 'source_files/combined.SLX_table_with_controls_031220.csv',
        BLOODSPOT = FALSE,
        COSMIC = TRUE
    )
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

# Filter rows that are on target.

filter.for.ontarget <- function(mutationTable)
{
    mutationTable %>% filter(ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)
}

##
# From TAPAS_functions.R
#
#' Annotate with locus error rate
#' Locus error rate = overall background error rate per locus, aggregated across control samples
annotate_with_locus_error_rate <- function(mutationTable,
                                           layoutTable,
                                           on_target,
                                           proportion_of_controls = 0.1,
                                           max_background_mean_AF = 0.01,
                                           is.blood_spot = FALSE,
                                           errorRateFile = NULL)
{
    if (is.blood_spot)
    {
        message("sWGS/blood spot mutationTable, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
        max_background_mean_AF <- 1
    }

    errorRateTable <- NULL

    if (!is.null(errorRateFile) && file.exists(errorRateFile))
    {
        # To speed up testing.
        warning("Reading error rate table from ", errorRateFile)
        errorRateTable <- readRDS(errorRateFile)
    }

    if (is.null(errorRateTable))
    {
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

        if (!is.null(errorRateFile))
        {
            saveRDS(errorRateTable, errorRateFile)
        }
    }

    errorRateTable <- errorRateTable %>%
        mutate(PROPORTION = N_SAMPLES_WITH_SIGNAL / N_SAMPLES,
               LOCUS_NOISE.PASS = PROPORTION < proportion_of_controls & BACKGROUND_AF < max_background_mean_AF,
               HAS_AF = BACKGROUND_AF > 0,
               LOCUS_NOISE.FAIL = !LOCUS_NOISE.PASS)

    counts <- errorRateTable %>%
        summarize(NON_ZERO_LOCI_PERCENT = sum(HAS_AF) / n() * 100,
                  LOCUS_NOISE_FAIL_PERCENT = sum(LOCUS_NOISE.FAIL) / n() * 100)

    ## TODO Copy and update plot code.
    ## That said, seems that's for a special case that doesn't apply for most work.

    proportionPositions <- errorRateTable %>%
        filter(PROPORTION < proportion_of_controls)

    backgroundAFPositions <- errorRateTable %>%
        filter(BACKGROUND_AF < max_background_mean_AF)

    extendedMutationTable <- mutationTable %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% proportionPositions$UNIQUE_POS & UNIQUE_POS %in% backgroundAFPositions$UNIQUE_POS)

    extendedMutationTable
}

##
# This function is built from the code after loading in INVAR3.R
#
# take only off_target bases, exclude INDELs
take.offtarget <- function(mutationTable, layoutTable, use_cosmic, is.blood_spot = FALSE)
{
    # Common function. Takes in a mutation table.
    # Adds MUT_SUM, DP_SUM, BACKGROUND_AF columns on rows grouped by
    # POOL, BARCODE, REF, ALT, TRINUCLEOTIDE
    groupAndSummarize <- function(mt)
    {
        mt %>%
            group_by(POOL, BARCODE, REF, ALT, TRINUCLEOTIDE) %>%
            summarize(MUT_SUM = sum(ALT_F) + sum(ALT_R),
                      DP_SUM = sum(DP),
                      BACKGROUND_AF= (sum(ALT_F) + sum(ALT_R)) / sum(DP),
                      .groups = 'drop')

    }

    mutationTable.off_target <- mutationTable %>%
        filter(!ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)

    if (!use_cosmic)
    {
        mutationTable.off_target <- mutationTable.off_target %>%
            filter(!COSMIC)
    }

    # error rate per locus filters
    backgroundErrorTable <- mutationTable.off_target %>%
        groupAndSummarize()

    locusErrorRateFile <- str_c('locus_error_rates.off_target.', ifelse(use_cosmic, 'cosmic', 'no_cosmic'), '.rds')

    mutations.off_target <- annotate_with_locus_error_rate(mutationTable.off_target,
                                                           layoutTable,
                                                           is.blood_spot = is.blood_spot,
                                                           on_target = FALSE,
                                                           errorRateFile = locusErrorRateFile)

    # annotate with both strand information
    mutations.off_target <- mutations.off_target %>%
        mutate(BOTH_STRANDS = (ALT_R > 0 & ALT_F > 0) | AF == 0)

    ## Calculate error rate with differnt settings

    # locus noise pass

    background_error.locus_noise <- mutations.off_target %>%
        filter(LOCUS_NOISE.PASS) %>%
        groupAndSummarize()

    # both strands

    background_error.both_strands <- mutations.off_target %>%
        filter(BOTH_STRANDS) %>%
        groupAndSummarize()

    # locus noise AND both strands

    background_error.locus_noise.both_strands <- mutations.off_target %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS) %>%
        groupAndSummarize()

    ## Assemble all error rates.

    error_rates <- list(pre_filter = backgroundErrorTable,
                        locus_noise_filter_only = background_error.locus_noise,
                        both_strands_only = background_error.both_strands,
                        locus_noise.both_strands = background_error.locus_noise.both_strands)

    error_rates
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    tumourMutationTable <- load.tumour.mutations.table(scriptArgs$TUMOUR_MUTATIONS_FILE)

    layoutTable <- load.layout.file(scriptArgs$LAYOUT_FILE)

    mutationTable <- readRDS(scriptArgs$MUTATIONS_TABLE_FILE)

    saveRDS(on_target, str_c(scriptArgs$TAPAS_SETTING, '.on_target.rds'))
    write_tsv(on_target, str_c(scriptArgs$TAPAS_SETTING, '.on_target.tsv'))

    off_target.cosmic <- take.offtarget(mutationTable, layoutTable, TRUE, scriptArgs$BLOODSPOT)
    saveRDS(off_target.cosmic, str_c(scriptArgs$TAPAS_SETTING, '.off_target.cosmic.error_rates.rds'))

    off_target.no_cosmic <- take.offtarget(mutationTable, layoutTable, FALSE, scriptArgs$BLOODSPOT)
    saveRDS(off_target.no_cosmic, str_c(scriptArgs$TAPAS_SETTING, '.off_target.no_cosmic.error_rates.rds'))
}

# Launch it.

invisible(main(parseOptions()))
