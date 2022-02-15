suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
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

    patientSummaryTable %>%
        write_csv("tumour_mutation_per_patient.csv")

    errorRatesINV042 <-
        calculateErrorRatesINV042(offTargetErrorRatesList[['PREFILTER']], layoutTable)

    calculateErrorRateSummary(errorRatesINV042) %>%
        write_csv("error_rate_summary.csv")

    sizeCharacterisationSummary <-
        calculateSizeCharacterisationSummary(sizeCharacterisationTable, layoutTable,
                                             study = scriptArgs$STUDY, roundTo = 5L)

    ifPatientData <-
        getIFPatientData(invarScoresTable, layoutTable, patientSummaryTable)

    exportCSV(ifPatientData$PATIENT_SPECIFIC_GLRT, 'patient_specific_GLRT.csv')
    exportCSV(ifPatientData$IF_PATIENT_DATA, 'IF_patient_data.csv')
    exportCSV(ifPatientData$THRESHOLD_EFFECTS, 'IR_threshold_effects.csv')

    ## Creating and saving plots.

    # plotting 3bp context

    ggsave(plot = cohortMutationContextPlot(contextMutationsTable, study = scriptArgs$STUDY),
           filename = "p1_cohort_mut_context.pdf",
           width = 6, height = 5)

    # plot the AF by mutation class

    ggsave(plot = cohortMutationAFByClassPlot(contextMutationsTable, study = scriptArgs$STUDY),
           filename = "p2_cohort_mut_AF_by_class.pdf",
           width = 6, height = 7)

    # plot mutations per patient captured and passing pipeline filters

    ggsave(plot = mutationsPerPatientPlot(patientSummaryTable, study = scriptArgs$STUDY),
           filename = "p3_cohort_mut_count_tumour.pdf",
           width = 6, height = 7)

    # plot mutation class distribution by cohort

    ggsave(plot = mutationClassByCohortPlot(contextMutationsTable, study = scriptArgs$STUDY),
           filename = "p4_cohort_mut_class.pdf",
           width = 6, height = 5)

    ## Case vs control error rates.

    ggsave(plot = backgroundErrorCaseControlPlot(errorRatesINV042, study = scriptArgs$STUDY),
           filename = "p5_background_error_comparison.case_vs_control.pdf",
           width = 6, height = 4)

    ## Summary plots

    ggsave(plot = backgroundErrorRatesPlot(errorRatesTable, study = scriptArgs$STUDY, familySize = scriptArgs$FAMILY_SIZE),
           filename = "p7_background_error_rates.pdf",
           width = 8, height = 4)

    ggsave(plot = errorRatePolishingComparisonPlot(errorRatesTable, "both_strands"),
           filename = "p9_error_rate_comparison.both_strands.pdf",
           width = 6, height = 3)

    ggsave(plot = errorRatePolishingComparisonPlot(errorRatesTable, "locus_noise"),
           filename = "p9_error_rate_comparison.locus_noise.pdf",
           width = 6, height = 3)

    ggsave(plot = errorRatePolishingComparisonPlot(errorRatesTable, "locus_noise.both_strands"),
           filename = "p9_error_rate_comparison.locus_noise_both_strands.pdf",
           width = 6, height = 3)

    ggsave(plot = filteringComparisonPlot(errorRatesTable, study = scriptArgs$STUDY, familySize = scriptArgs$FAMILY_SIZE),
           filename = "p20_filtering_comparison.pdf",
           width = 6, height = 4)

    ## Outlier suppression plots

    ggsave(plot = summaryCohortPlot(mutationsTable, study = scriptArgs$STUDY),
           filename = "p6_os_data_retained.pdf",
           width = 5, height = 6)

    ggsave(plot = backgroundPolishingPlot(mutationsTable, study = scriptArgs$STUDY, errorSuppression = scriptArgs$ERROR_SUPPRESSION, outlierSuppression = scriptArgs$OUTLIER_SUPPRESSION),
           filename = "p8_background_polishing.pdf",
           width = 15, height = 12)

    ## Tumour AF in observed and unobserved loci.

    ggsave(plot = tumourAFInLociPlot(mutationsTable, study = scriptArgs$STUDY),
           filename = "p10_tumour_AF_for_observed_non_observed_loci.pdf",
           width = 5, height = 4)

    ## Fragment size per cohort, with different levels of error-suppression

    ggsave(plot = fragmentSizePerCohortPlot(sizeCharacterisationSummary, study = scriptArgs$STUDY),
           filename = "p11_size_comparison.pdf",
           width = 6, height = 5)

    ## Enrichment level

    ggsave(plot = enrichmentLevelPlot(sizeCharacterisationSummary, study = scriptArgs$STUDY),
           filename = "p12_enrichment_ratios.pdf",
           width = 6, height = 5)

    ## Receiver Operating Characteristic Plots

    ggsave(plot = receiverOperatingCharacteristicPlot(invarScoresTable, layoutTable, TRUE,
                                                      study = scriptArgs$STUDY,
                                                      familySize = scriptArgs$FAMILY_SIZE),
           filename = "p13a_receiver_operating_characteristic.pdf",
           width = 4, height = 3)

    ggsave(plot = receiverOperatingCharacteristicPlot(invarScoresTable, layoutTable, FALSE,
                                                      study = scriptArgs$STUDY,
                                                      familySize = scriptArgs$FAMILY_SIZE),
           filename = "p13b_receiver_operating_characteristic.no_size.pdf",
           width = 4, height = 3)

    ## IR (depth) to IMAF plot

    ggsave(plot = depthToIMAFPlot(ifPatientData$IF_PATIENT_DATA),
           filename = "p14_IR_vs_IMAF.pdf",
           width = 6, height = 4)

    ## Waterfall plot with detectable vs non detectable using dPCR

    ggsave(plot = detectableWaterfallPlot(ifPatientData$PATIENT_SPECIFIC_GLRT,
                                          layoutTable,
                                          patientSummaryTable,
                                          study = scriptArgs$STUDY),
           filename = "p15_waterfall_IMAF.pdf",
           width = 10, height = 5)

    ## Waterfall plot of cancer genomes

    ggsave(plot = cancerGenomesWaterfallPlot(ifPatientData$PATIENT_SPECIFIC_GLRT,
                                             layoutTable,
                                             patientSummaryTable,
                                             study = scriptArgs$STUDY),
           filename = "p16_waterfall_cancer_genomes.pdf",
           width = 10, height = 5)
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR/analysis')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
