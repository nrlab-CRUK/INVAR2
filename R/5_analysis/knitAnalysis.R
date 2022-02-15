suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plotROC))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))
source(str_c(Sys.getenv('INVAR_HOME'), '/R/5_analysis/analysisCalculations.R'))
source(str_c(Sys.getenv('INVAR_HOME'), '/R/5_analysis/analysisPlots.R'))


##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"

    options_list <- list(
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) with outlier suppression flags",
                    default=defaultMarker),
        make_option(c("--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("--error-rates"), type="character", metavar="file",
                    dest="ERROR_RATES_FILE", help="The on target background error rates.",
                    default=defaultMarker),
        make_option(c("--off-target-error-rates"), type="character", metavar="file",
                    dest="OFF_TARGET_ERROR_RATES_FILE", help="The off target error rates.",
                    default=defaultMarker),
        make_option(c("--size-characterisation"), type="character", metavar="file",
                    dest="SIZE_CHARACTERISATION_FILE", help="The size characterisation file.",
                    default=defaultMarker),
        make_option(c("--invar-scores"), type="character", metavar="file",
                    dest="INVAR_SCORES_FILE", help="The INVAR scores file.",
                    default=defaultMarker),
        make_option(c("--study"), type="character", metavar="string",
                    dest="STUDY", help="The study name",
                    default=defaultMarker),
        make_option(c("--error-suppression"), type="character", metavar="string",
                    dest="ERROR_SUPPRESSION", help="The error suppression string",
                    default=defaultMarker),
        make_option(c("--family-size"), type="integer", metavar="integer",
                    dest="FAMILY_SIZE", help="The family size",
                    default=defaultMarker),
        make_option(c("--outlier-suppression"), type="double", metavar="number",
                    dest="OUTLIER_SUPPRESSION", help="The outlier suppression threshold",
                    default=0.05))

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
    testhome <- str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/')
    base <- str_c(testhome, 'analysis/source/')

    list(
        MUTATIONS_TABLE_FILE = str_c(base, 'mutation_table.outliersuppressed.rds'),
        ERROR_RATES_FILE = str_c(base, 'background_error_rates.rds'),
        OFF_TARGET_ERROR_RATES_FILE = str_c(base, 'error_rates.off_target.no_cosmic.rds'),
        SIZE_CHARACTERISATION_FILE = str_c(base, 'size_characterisation.rds'),
        INVAR_SCORES_FILE = str_c(base, 'invar_scores.rds'),
        LAYOUT_FILE = str_c(testhome, 'invar_source/combined.SLX_table_with_controls_031220.v2.csv'),
        STUDY = 'PARADIGM',
        ERROR_SUPPRESSION = 'f0.9_s2',
        FAMILY_SIZE = 2L,
        OUTLIER_SUPPRESSION = 0.05
    )
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    layoutTable <-
        loadLayoutTable(scriptArgs$LAYOUT_FILE)

    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    errorRatesTable <- readRDS(scriptArgs$ERROR_RATES_FILE) %>%
        mutate(MUTATION_CLASS = str_c(REF, ALT, sep = '/'))

    offTargetErrorRatesList <- readRDS(scriptArgs$OFF_TARGET_ERROR_RATES_FILE)

    sizeCharacterisationTable <- readRDS(scriptArgs$SIZE_CHARACTERISATION_FILE)

    invarScoresTable <- readRDS(scriptArgs$INVAR_SCORES_FILE) %>%
        mutate(PATIENT_SPECIFIC = PATIENT == PATIENT_MUTATION_BELONGS_TO,
               POOL_BARCODE = str_c(POOL, BARCODE, sep = "_"))

    assert_that(length(offTargetErrorRatesList) == 4 && all(names(offTargetErrorRatesList) %in% c('PREFILTER', 'LOCUS_NOISE', 'BOTH_STRANDS', 'LOCUS_NOISE.BOTH_STRANDS')),
                msg = "Off target error rates list does not contain the tables expected.")


    # Manipulation and further calculations.

    contextMutationsTable <- mutationsTable %>%
        filter(PATIENT_SPECIFIC & LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & OUTLIER.PASS) %>%
        group_by(PATIENT, UNIQUE_POS) %>%
        slice_head(n = 1) %>%
        ungroup()

    patientSummaryTable <- contextMutationsTable %>%
        group_by(PATIENT) %>%
        summarise(MUTATIONS = n_distinct(PATIENT, UNIQUE_POS), .groups = "drop") %>%
        arrange(PATIENT, MUTATIONS)

    errorRatesINV042 <-
        calculateErrorRatesINV042(offTargetErrorRatesList[['PREFILTER']], layoutTable)

    sizeCharacterisationSummary <-
        calculateSizeCharacterisationSummary(sizeCharacterisationTable, layoutTable,
                                             study = scriptArgs$STUDY, roundTo = 5L)

    ifPatientData <-
        getIFPatientData(invarScoresTable, layoutTable, patientSummaryTable)

    annotatedPatientSpecificGLRT <-
        annotatePatientSpecificGLRT(ifPatientData$PATIENT_SPECIFIC_GLRT, layoutTable, patientSummaryTable)

    ## Create the report.

    rmarkdown::render(input = str_c(Sys.getenv('INVAR_HOME'), '/R/5_analysis/analysis.Rmd'),
                      knit_root_dir = str_c(getwd(), '/knitting'),
                      intermediates_dir = str_c(getwd(), '/intermediates'),
                      output_format = rmarkdown::html_document(),
                      output_dir = getwd(),
                      output_file = str_c(study, "_invar2_analysis.html")
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR/analysis')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
