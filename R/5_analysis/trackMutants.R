suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))
source(str_c(Sys.getenv('INVAR_HOME'), '/R/5_analysis/analysisCalculations.R'))


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
        make_option(c("--invar-scores"), type="character", metavar="file",
                    dest="INVAR_SCORES_FILE", help="The INVAR scores file.",
                    default=defaultMarker),
        make_option(c("--tumour-mutations"), type="character", metavar="file",
                    dest="TUMOUR_MUTATIONS_FILE", help="The source patient mutations file",
                    default=defaultMarker),
        make_option(c("--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("--study"), type="character", metavar="string",
                    dest="STUDY", help="The study name",
                    default=defaultMarker))

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
    base <- str_c(testhome, 'trackMutants/source/')

    list(
        MUTATIONS_TABLE_FILE = str_c(base, 'mutation_table.outliersuppressed.rds'),
        INVAR_SCORES_FILE = str_c(base, 'invar_scores.rds'),
        TUMOUR_MUTATIONS_FILE = str_c(testhome, 'invar_source/PARADIGM_mutation_list_full_cohort_hg19.v2.csv'),
        LAYOUT_FILE = str_c(testhome, 'invar_source/combined.SLX_table_with_controls_031220.v2.csv'),
        STUDY = 'PARADIGM'
    )
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    layoutTable <-
        loadLayoutTable(scriptArgs$LAYOUT_FILE)

    layoutTableTimepoint <- layoutTable %>%
        select(POOL, BARCODE, TIMEPOINT)

    tumourMutationTableSummary <-
        loadTumourMutationsTable(scriptArgs$TUMOUR_MUTATIONS_FILE) %>%
        group_by(PATIENT) %>%
        summarise(INITIAL_MUTATIONS = n(), .groups = "drop")

    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns() %>%
        mutate(UNIQUE_PATIENT_POS = str_c(UNIQUE_POS, UNIQUE_ALT, sep = '_'))

    # Only interested in patient specific rows from INVAR scores.
    # LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS are all true for the current
    # set up, so they don't need to be included in this table.

    invarScoresTable <- readRDS(scriptArgs$INVAR_SCORES_FILE) %>%
        mutate(PATIENT_SPECIFIC = PATIENT == PATIENT_MUTATION_BELONGS_TO) %>%
        adjustInvarScores(layoutTable) %>%
        filter(PATIENT_SPECIFIC) %>%
        select(POOL, BARCODE, PATIENT, USING_SIZE, OUTLIER.PASS, DETECTION)

    mutationTracking <- mutationsTable %>%
        mutate(UNIQUE_IF_SPECIFIC_MUTANT = ifelse(PATIENT_SPECIFIC & MUTANT, UNIQUE_PATIENT_POS, NA),
               UNIQUE_IF_MUTANT = ifelse(MUTANT, UNIQUE_PATIENT_POS, NA),
               UNIQUE_IF_NOT_MUTANT = ifelse(!MUTANT, UNIQUE_PATIENT_POS, NA)) %>%
        group_by(POOL, BARCODE, PATIENT, CASE_OR_CONTROL) %>%
        summarise(LOCUS_NOISE.PASS = sum(LOCUS_NOISE.PASS & PATIENT_SPECIFIC),
                  PATIENT_SPECIFIC_WITH_MUTATION = n_distinct(UNIQUE_IF_SPECIFIC_MUTANT, na.rm = TRUE),
                  MUTANTS_PATIENT_SPECIFIC = n_distinct(UNIQUE_IF_MUTANT, na.rm = TRUE),
                  MUTANTS_NON_SPECIFIC = n_distinct(UNIQUE_IF_NOT_MUTANT, na.rm = TRUE),
                  SPECIFIC_OUTLIER.PASS = sum(PATIENT_SPECIFIC & OUTLIER.PASS),
                  SPECIFIC.PASS = sum(PATIENT_SPECIFIC & LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS & OUTLIER.PASS & MUTATION_SUM > 0),
                  NON_SPECIFIC.PASS = sum(!PATIENT_SPECIFIC & LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS & OUTLIER.PASS & MUTATION_SUM > 0),
                  .groups = "drop") %>%
        left_join(layoutTableTimepoint, by = c("POOL", "BARCODE")) %>%
        left_join(tumourMutationTableSummary, by = "PATIENT") %>%
        full_join(invarScoresTable, by = c('POOL', 'BARCODE', 'PATIENT')) %>%
        # mutate(INITIAL_MUTATIONS = ifelse(CASE_OR_CONTROL == 'case', INITIAL_MUTATIONS, 0)) %>%
        mutate(USING_SIZE = ifelse(USING_SIZE, 'WITH_SIZE', 'NO_SIZE'),
               OUTLIER.PASS = ifelse(OUTLIER.PASS, 'PASS', 'FAIL')) %>%
        pivot_wider(names_from = c(USING_SIZE, OUTLIER.PASS),
                    names_glue = "{USING_SIZE}.OUTLIER_{OUTLIER.PASS}",
                    values_from = DETECTION) %>%
        select(POOL, BARCODE, PATIENT, TIMEPOINT, CASE_OR_CONTROL,
               INITIAL_MUTATIONS, LOCUS_NOISE.PASS,
               PATIENT_SPECIFIC_WITH_MUTATION,
               MUTANTS_PATIENT_SPECIFIC, MUTANTS_NON_SPECIFIC,
               SPECIFIC_OUTLIER.PASS, SPECIFIC.PASS, NON_SPECIFIC.PASS,
               WITH_SIZE.OUTLIER_PASS, NO_SIZE.OUTLIER_PASS,
               WITH_SIZE.OUTLIER_FAIL, NO_SIZE.OUTLIER_FAIL) %>%
        arrange(POOL, BARCODE, PATIENT, TIMEPOINT, CASE_OR_CONTROL)

    mutationTracking %>%
        arrangeMutationTableForExport() %>%
        exportCSV('mutations_tracking.csv')

    mutationTracking
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    summary <- invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
