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

    options_list <- list()

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
        stop("ERROR: Need at least one file from generalisedLikelihoodRatioTest for combination.")
    }

    scriptOptions <- list(GLRT_TABLE_FILES = opts$args)
    for (v in names(opts$options))
    {
        scriptOptions[v] = opts$options[v]
    }

    scriptOptions
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    glrtResultsList <- lapply(scriptArgs$GLRT_TABLE_FILES, readRDS)

    glrtTable <-
        bind_rows(glrtResultsList) %>%
        arrange(SAMPLE_ID, PATIENT_MUTATION_BELONGS_TO,
                LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS, ITERATION, USING_SIZE)

    saveRDS(glrtTable, 'invar_scores.rds')
}

# Launch it.

invisible(main(parseOptions()))
