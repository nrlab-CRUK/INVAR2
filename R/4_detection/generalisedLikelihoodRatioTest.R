suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
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
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by outlierSuppression.R",
                    default=defaultMarker),
        make_option(c("--size-characterisation"), type="character", metavar="file",
                    dest="SIZE_CHARACTERISATION_FILE", help="The size characterisation summary file (RDS)",
                    default=defaultMarker),
        make_option(c("--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."))

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
    base <- str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/generalisedLikelihoodRatioTest/source')

    list(
        MUTATIONS_TABLE_FILE = str_c(base, '/SLX-19721_SXTLI001.os.rds'),
        SIZE_CHARACTERISATION_FILE = str_c(base, '/size_characterisation.rds'),
        BLOODSPOT = FALSE
    )
}


##
# From TAPAS_functions.R, originally calculate_IMAFv2.
#
calculateIMAFv2 <- function(mutationsTable, bloodspot)
{
    assert_that(is.logical(bloodspot), msg = "bloodspot argument must be a logical")

    mutationsTable.flat <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS, BOTH_STRANDS) %>%
        select(-SIZE, -MUTANT) %>%
        distinct()

    # When it's blood spot data, no need to use outlier suppression as depth is either 1 or 0
    # When not blood spot, also filter on outlier suppression (not an outlier).
    if (!bloodspot)
    {
        mutationsTable.flat <- mutationsTable.flat %>%
            filter(OUTLIER.PASS)
    }

    summary <- mutationsTable.flat %>%
        group_by(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO, MUTATION_CLASS, TRINUCLEOTIDE) %>%
        summarise(TOTAL_DP = sum(REF_F + REF_R + ALT_F + ALT_R),
                  MUTATION_SUM = sum(MUTATION_SUM),
                  MEAN_AF = weighted.mean(AF, DP),
                  BACKGROUND_AF.TRINUCLEOTIDE = first(BACKGROUND_AF),
                  .groups = "drop") %>%
        group_by(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO) %>%
        mutate(MEAN_AF.BS_TRINUCLEOTIDE = pmax(0, MEAN_AF - BACKGROUND_AF.TRINUCLEOTIDE)) %>%
        summarise(IMAFV2 = signif(weighted.mean(MEAN_AF.BS_TRINUCLEOTIDE, TOTAL_DP), 4),
                  .groups = "drop")

    summary
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    # TESTING
    #mutationTable.one <- mutationTable %>% filter(PATIENT_MUTATION_BELONGS_TO == 'PARA_002')

    #sizeCharacterisationTable <- readRDS(scriptArgs$SIZE_CHARACTERISATION_FILE)

    IMAFv2 <- calculateIMAFv2(mutationTable, bloodspot = scriptArgs$BLOODSPOT)

    IMAFv2
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
