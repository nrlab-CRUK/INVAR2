# This is basically annotate_with_outlier_suppression.R in the original pipeline

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
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by part one of the pipeline",
                    default=defaultMarker),
        make_option(c("--threads"), type="integer", metavar="integer",
                    dest="THREADS", help="The number of cores to use to process the input files.",
                    default=1))

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
        stop("ERROR: Need at least one file from outlierSuppression.R.")
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
    base <- str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/annotateMutationsWithOutlierSuppression/source')

    files = list.files(base, pattern = "\\.os.rds$", full.names = TRUE)

    list(
        MUTATIONS_TABLE_FILE = str_c(base, '/mutation_table.on_target.rds'),
        MUTATIONS_TABLE_FILES = files,
        THREADS = 4
    )
}


##
# Extract just the outlier pass flag from the mutation table file from
# outlierSuppression.R.
#

getMutationOutlierFlags <- function(mutationTableFileName)
{
    message("Filtering ", mutationTableFileName)

    # This creates the equivalent of the "nonzero" files made by the
    # original outlier_suppression.R script. It's unique positions
    # in each source file that are mutants, selecting an example.

    mutationTable <- readRDS(mutationTableFileName)

    reducedTable <- mutationTable %>%
        filter(MUTANT) %>%
        group_by(POOL, BARCODE, CHROM, POS) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        select(POOL, BARCODE, CHROM, POS, OUTLIER.PASS)

    reducedTable
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE)

    outlierFlagTables <-
        mclapply(scriptArgs$MUTATIONS_TABLE_FILES, getMutationOutlierFlags, mc.cores = scriptArgs$THREADS)

    outlierFlagTable <-
        bind_rows(outlierFlagTables) %>%
        arrange(POOL, BARCODE, CHROM, POS)

    saveRDS(outlierFlagTable, "outlierflags.rds")

    outlierFlagTable %>%
        filter(!OUTLIER.PASS) %>%
        exportTSV("outlierflags.fail.tsv")

    outlierFlagTable %>%
        filter(OUTLIER.PASS) %>%
        exportTSV("outlierflags.pass.tsv")

    # Join tables with the outlier flag tables. Where a row doesn't match it will
    # have NA as the value. These can be changed to TRUE.
    # I suppose those rows that have no match are not found to be outliers, so the
    # default of passing is sensible.

    mutationTable.withOutlier <- mutationTable %>%
        left_join(outlierFlagTable, by = c('POOL', 'BARCODE', 'CHROM', 'POS')) %>%
        mutate(OUTLIER.PASS = ifelse(is.na(OUTLIER.PASS), TRUE, OUTLIER.PASS))

    mutationTable.withOutlier %>%
        arrangeMutationTableForExport() %>%
        saveRDS('mutation_table.with_outliers.rds')
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
