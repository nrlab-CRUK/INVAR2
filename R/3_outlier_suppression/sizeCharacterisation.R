# This is basically size_annotation3.R in the original pipeline

suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
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
        make_option(c("--threads"), type="integer", metavar="integer",
                    dest="THREADS", help="The number of cores to use to process the input files.",
                    default=1L))

    opts <- OptionParser(option_list=options_list, usage="%prog [options]") %>%
        parse_args(positional_arguments = TRUE)

    missingRequired <- which(opts$options == defaultMarker)
    if (any(missingRequired))
    {
        missing <- names(opts$options)[missingRequired]
        stop(str_c("ERROR: One or more required arguments have not been provided: ", str_c(missing, collapse = ' ')))
    }

    if (length(opts$args) < 1)
    {
        stop("ERROR: Need at least one file from outlierSuppression.R to get size characteristics from.")
    }

    scriptOptions <- list(MUTATIONS_TABLE_FILES = opts$args)
    for (v in names(opts$options))
    {
        scriptOptions[v] = opts$options[v]
    }

    scriptOptions
}

# Test options for my (Rich) local set up in RStudio.

richTestOptions <- function()
{
    files = list.files(str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/sizeCharacterisation/source'),
                       pattern = "\\.rds$", full.names = TRUE)

    list(
        MUTATIONS_TABLE_FILES = files,
        THREADS = 1L
    )
}


##
# Calculate the size characterisation per file.
# This was the body of sizeAnnotation3.R
#

calculateSizeCharacteristics <- function(mutationTableFileName)
{
    message("Summarising ", mutationTableFileName)

    mutationTable <- readRDS(mutationTableFileName) %>%
        mutate(PATIENT_SPECIFIC = PATIENT == PATIENT_MUTATION_BELONGS_TO)

    sizeCharacteristicsTable <- mutationTable %>%
        filter(OUTLIER.PASS & BOTH_STRANDS.PASS & LOCUS_NOISE.PASS) %>%
        group_by(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO, CASE_OR_CONTROL, PATIENT_SPECIFIC, MUTANT, SIZE) %>%
        summarise(TOTAL = n(), .groups = "drop")

    sizeCharacteristicsTable
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    sizeTables <-
        mclapply(scriptArgs$MUTATIONS_TABLE_FILES, calculateSizeCharacteristics, mc.cores = scriptArgs$THREADS)

    sizeCharacteristicsTable <-
        bind_rows(sizeTables) %>%
        arrange(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO, MUTANT, SIZE)

    sizeCharacteristicsSummary <- sizeCharacteristicsTable %>%
        group_by(PATIENT_SPECIFIC, CASE_OR_CONTROL, MUTANT, SIZE) %>%
        summarise(COUNT = sum(TOTAL), .groups = "drop") %>%
        group_by(PATIENT_SPECIFIC, CASE_OR_CONTROL, MUTANT) %>%
        mutate(TOTAL = sum(COUNT)) %>%
        ungroup() %>%
        mutate(PROPORTION = COUNT / TOTAL) %>%
        arrange(PATIENT_SPECIFIC, CASE_OR_CONTROL, MUTANT, SIZE)

    sizeCharacteristicsTable %>%
        saveRDS("size_characterisation.all.rds")

    sizeCharacteristicsSummary %>%
        saveRDS("size_characterisation.rds")
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
