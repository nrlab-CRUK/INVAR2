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
    make_option(c("--tumour-mutations"), type="character", metavar="file",
                dest="TUMOUR_MUTATIONS_FILE", help="The source patient mutations file",
                default=defaultMarker),
    make_option(c("--layout"), type="character", metavar="file",
                dest="LAYOUT_FILE", help="The sequencing layout file",
                default=defaultMarker),
    make_option(c("--error-rates"), type="character", metavar="file",
                dest="ERROR_RATES_FILE", help="The on target background error rates.",
                default=defaultMarker),
    make_option(c("--on-target-locus-error-rates"), type="character", metavar="file",
                dest="ON_TARGET_LOCUS_ERROR_RATES_FILE", help="The on target locus error rates.",
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
    make_option(c("--tapas"), type="character", metavar="string",
                dest="TAPAS_SETTING", help="The TAPAS setting",
                default=defaultMarker),
    make_option(c("--error-suppression"), type="character", metavar="string",
                dest="ERROR_SUPPRESSION", help="The error suppression string",
                default=defaultMarker),
    make_option(c("--family-size"), type="integer", metavar="integer",
                dest="FAMILY_SIZE", help="The family size",
                default=defaultMarker),
    make_option(c("--outlier-suppression"), type="double", metavar="number",
                dest="OUTLIER_SUPPRESSION", help="The outlier suppression threshold",
                default=0.05),
    make_option(c("--max-background-allele-frequency"), type="double", metavar="num",
                dest="MAX_BACKGROUND_ALLELE_FREQUENCY", help="Filter loci with a background allele frequency in controls greater than this value",
                default=0.01),
    make_option(c("--allele-frequency-threshold"), type="double", metavar="num",
                dest="ALLELE_FREQUENCY_THRESHOLD", help="Maximum allele frequency value for acceptable samples.",
                default=0.01),
    make_option(c("--max-mutant-reads"), type="integer", metavar="integer",
                dest="MAXIMUM_MUTANT_READS", help="Maximum number of reads acceptable for outlier suppression trying to detect MRD",
                default=10),
    make_option(c("--min-N-IR"), type="integer", metavar="integer",
                dest="MINIMUM_N_INFORMATIVE_READS", help="Minimum number of informative reads per sample for it to be considered as part of the cohort",
                default=10),
    make_option(c("--score-specificity"), type="double", metavar="number",
                dest="SCORE_SPECIFICITY", help="Score specificity for ROC plot.",
                default=0.95))
  
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


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
  assert_that(file.exists(scriptArgs$LAYOUT_FILE), msg = str_c(scriptArgs$LAYOUT_FILE, " does not exist."))
  assert_that(file.exists(scriptArgs$TUMOUR_MUTATIONS_FILE), msg = str_c(scriptArgs$TUMOUR_MUTATIONS_FILE, " does not exist."))
  assert_that(file.exists(scriptArgs$MUTATIONS_TABLE_FILE), msg = str_c(scriptArgs$MUTATIONS_TABLE_FILE, " does not exist."))
  assert_that(file.exists(scriptArgs$ERROR_RATES_FILE), msg = str_c(scriptArgs$ERROR_RATES_FILE, " does not exist."))
  assert_that(file.exists(scriptArgs$OFF_TARGET_ERROR_RATES_FILE), msg = str_c(scriptArgs$OFF_TARGET_ERROR_RATES_FILE, " does not exist."))
  assert_that(file.exists(scriptArgs$SIZE_CHARACTERISATION_FILE), msg = str_c(scriptArgs$SIZE_CHARACTERISATION_FILE, " does not exist."))
  assert_that(file.exists(scriptArgs$INVAR_SCORES_FILE), msg = str_c(scriptArgs$INVAR_SCORES_FILE, " does not exist."))
  
  layoutTable <-
    loadLayoutTable(scriptArgs$LAYOUT_FILE)
  
  tumourMutationsTable <-
    loadTumourMutationsTable(scriptArgs$TUMOUR_MUTATIONS_FILE)
  
  mutationsTable <-
    readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
    addMutationTableDerivedColumns()
  
  errorRatesTable <- readRDS(scriptArgs$ERROR_RATES_FILE) %>%
    mutate(MUTATION_CLASS = str_c(REF, ALT, sep = '/'))
  
  onTargetLocusErrorRatesTable <- readRDS(scriptArgs$ON_TARGET_LOCUS_ERROR_RATES_FILE)
  
  offTargetErrorRatesList <- readRDS(scriptArgs$OFF_TARGET_ERROR_RATES_FILE)
  
  sizeCharacterisationTable <- readRDS(scriptArgs$SIZE_CHARACTERISATION_FILE)
  
  invarScoresTable <- readRDS(scriptArgs$INVAR_SCORES_FILE) %>%
    mutate(PATIENT_SPECIFIC = PATIENT == PATIENT_MUTATION_BELONGS_TO)
  
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
  
  inputMutationsTable <- mutationsTable %>%
    filter(PATIENT_SPECIFIC) %>%
    group_by(PATIENT) %>%
    summarise(INPUT_MUTATIONS = n_distinct(PATIENT, UNIQUE_POS), .groups = "drop") %>%
    arrange(PATIENT, INPUT_MUTATIONS)
  
  errorRatesINV042 <-
    calculateErrorRatesINV042(offTargetErrorRatesList[['PREFILTER']], layoutTable)
  
  sizeCharacterisationSummary <-
    calculateSizeCharacterisationSummary(sizeCharacterisationTable, layoutTable,
                                         study = scriptArgs$STUDY, roundTo = 5L)
  
  ifPatientData <-
    tryCatch(
      {
        getIFPatientData(invarScoresTable, layoutTable, patientSummaryTable, scriptArgs$SCORE_SPECIFICITY, scriptArgs$MINIMUM_N_INFORMATIVE_READS , scriptArgs$MAX_BACKGROUND_ALLELE_FREQUENCY)
      },
      error = function(cond)
      {
        warning(geterrmessage())
        NULL
      })
  
  annotatedPatientSpecificGLRT <- NULL
  if (!is.null(ifPatientData))
  {
    annotatedPatientSpecificGLRT <-
      annotatePatientSpecificGLRT(ifPatientData$PATIENT_SPECIFIC_GLRT, layoutTable, patientSummaryTable, scriptArgs$MINIMUM_N_INFORMATIVE_READS)
  }
  
  mutationTrackingTable <-
    mutationTracking(mutationsTable, layoutTable, tumourMutationsTable, invarScoresTable, scriptArgs$SCORE_SPECIFICITY)
  
  ## Saving some tables as part of analysis.
  
  exportCSV(patientSummaryTable, "tumour_mutation_per_patient.csv")
  
  exportCSV(calculateErrorRateSummary(errorRatesINV042), "error_rate_summary.csv")
  
  if (!is.null(ifPatientData))
  {
#    exportCSV(ifPatientData$PATIENT_SPECIFIC_GLRT, 'patient_specific_GLRT.csv')
#    exportCSV(ifPatientData$IF_PATIENT_DATA, 'IF_patient_data.csv')
    exportCSV(ifPatientData$THRESHOLD_EFFECTS, 'IR_threshold_effects.csv')
  }
  
  if (!is.null(annotatedPatientSpecificGLRT))
  {
    annotatedPatientSpecificGLRT <- annotatedPatientSpecificGLRT %>%
      mutate(UNIQUE_MOLECULES = DP / MUTATIONS,
             NG_ON_SEQ = UNIQUE_MOLECULES / 300,
             LOW_SENSITIVITY = DP < scriptArgs$MINIMUM_N_INFORMATIVE_READS & !DETECTED.WITH_SIZE)
    exportCSV(annotatedPatientSpecificGLRT, 'Results_summary.csv')
  }
  
  mutationTrackingTable %>%
    arrangeMutationTableForExport() %>%
    exportCSV('mutations_tracking.csv')
  
  ## Creating plots.
  
  plots <- list()
  
  # On target locus error rates.
  plots$P0 <- onTargetLocusErrorRatePlot(onTargetLocusErrorRatesTable,
                                         study = scriptArgs$STUDY,
                                         tapasSetting = scriptArgs$TAPAS_SETTING)
  
  # plotting 3bp context
  plots$P1 <- cohortMutationContextPlot(contextMutationsTable, study = scriptArgs$STUDY)
  
  # plot the AF by mutation class
  plots$P2 <- cohortMutationAFByClassPlot(contextMutationsTable, study = scriptArgs$STUDY)
  
  # plot mutations per patient captured and passing pipeline filters
  plots$P3 <- mutationsPerPatientPlot(patientSummaryTable, study = scriptArgs$STUDY)
  
  # plot mutations per patient captured and passing pipeline filters
  plots$P3_1 <- inUsedMutationsPerPatientPlot(patientSummaryTable, inputMutationsTable, study = scriptArgs$STUDY)
  
  # plot mutation class distribution by cohort
  plots$P4 <- mutationClassByCohortPlot(contextMutationsTable, study = scriptArgs$STUDY)
  
  # Case vs control error rates.
  plots$P5 <- backgroundErrorCaseControlPlot(errorRatesINV042, study = scriptArgs$STUDY)
  
  # Summary plots
  
  plots$P7 <- backgroundErrorRatesPlot(errorRatesTable,
                                       study = scriptArgs$STUDY,
                                       familySize = scriptArgs$FAMILY_SIZE)
  
  plots$P9a <- errorRatePolishingComparisonPlot(errorRatesTable, "both_strands")
  plots$P9b <- errorRatePolishingComparisonPlot(errorRatesTable, "locus_noise")
  plots$P9c <- errorRatePolishingComparisonPlot(errorRatesTable, "locus_noise.both_strands")
  
  plots$P20 <- filteringComparisonPlot(errorRatesTable,
                                       study = scriptArgs$STUDY,
                                       familySize = scriptArgs$FAMILY_SIZE)
  
  # Outlier suppression plots
  
  plots$P6 <- summaryCohortPlot(mutationsTable, study = scriptArgs$STUDY)
  
  plots$P8 <- backgroundPolishingPlot(mutationsTable, layoutTable,
                                      study = scriptArgs$STUDY,
                                      errorSuppression = scriptArgs$ERROR_SUPPRESSION,
                                      outlierSuppression = scriptArgs$OUTLIER_SUPPRESSION,
                                      allele_frequency_threshold = scriptArgs$ALLELE_FREQUENCY_THRESHOLD,
                                      MAXIMUM_MUTANT_READS = scriptArgs$MAXIMUM_MUTANT_READS)
  
  # Tumour AF in observed and unobserved loci.
  
  plots$P10 <- tumourAFInLociPlot(mutationsTable, study = scriptArgs$STUDY)
  
  # Fragment size per cohort, with different levels of error-suppression
  
  plots$P11 <- fragmentSizePerCohortPlot(sizeCharacterisationSummary, study = scriptArgs$STUDY)
  
  # Enrichment level
  
  plots$P12 <- enrichmentLevelPlot(sizeCharacterisationSummary, study = scriptArgs$STUDY)
  
  # Receiver Operating Characteristic Plots
  
  plots$P13a <- receiverOperatingCharacteristicPlot(invarScoresTable, layoutTable,
                                                    withSizes = TRUE,
                                                    study = scriptArgs$STUDY,
                                                    familySize = scriptArgs$FAMILY_SIZE,
                                                    scoreSpecificity = scriptArgs$SCORE_SPECIFICITY, 
                                                    MINIMUM_N_INFORMATIVE_READS = scriptArgs$MINIMUM_N_INFORMATIVE_READS, 
                                                    MAX_BACKGROUND_ALLELE_FREQUENCY = scriptArgs$MAX_BACKGROUND_ALLELE_FREQUENCY)
  
  plots$P13b <- receiverOperatingCharacteristicPlot(invarScoresTable, layoutTable,
                                                    withSizes = FALSE,
                                                    study = scriptArgs$STUDY,
                                                    familySize = scriptArgs$FAMILY_SIZE,
                                                    scoreSpecificity = scriptArgs$SCORE_SPECIFICITY,
                                                    MINIMUM_N_INFORMATIVE_READS = scriptArgs$MINIMUM_N_INFORMATIVE_READS, 
                                                    MAX_BACKGROUND_ALLELE_FREQUENCY = scriptArgs$MAX_BACKGROUND_ALLELE_FREQUENCY)
  
  # IR (depth) to IMAF plot
  
  plots$P14 <- NULL
  if (!is.null(ifPatientData))
  {
    plots$P14 <- depthToIMAFPlot(ifPatientData$IF_PATIENT_DATA)
  }
  
  plots$P15 <- NULL
  plots$P16 <- NULL
  plots$P17 <- NULL
  if (!is.null(annotatedPatientSpecificGLRT))
  {
    # Waterfall plot with detectable vs non detectable using dPCR
    
    plots$P15 <- detectableWaterfallPlot(annotatedPatientSpecificGLRT, study = scriptArgs$STUDY)
    
    # Waterfall plot of cancer genomes
    
    plots$P16 <- cancerGenomesWaterfallPlot(annotatedPatientSpecificGLRT, study = scriptArgs$STUDY)
    
    # dPCR comparision plot
    
    plots$P17 <- dpcrComparisionPlot(annotatedPatientSpecificGLRT, study = scriptArgs$STUDY)
  }
  
  
  ## Save the plots as individual files.
  
  plots$P0 <-
    savePlotSafely(plot = plots$P0,
                   filename = "p0_on_target_locus_error_rates.pdf",
                   width = 11, height = 7)
  
  plots$P1 <-
    savePlotSafely(plot = plots$P1,
                   filename = "p1_cohort_mut_context.pdf",
                   width = 6, height = 5)
  
  # plot the AF by mutation class
  
  plots$P2 <-
    savePlotSafely(plot = plots$P2,
                   filename = "p2_cohort_mut_AF_by_class.pdf",
                   width = 6, height = 7)
  
  # plot mutations per patient captured and passing pipeline filters
  
  plots$P3 <-
    savePlotSafely(plot = plots$P3,
                   filename = "p3_cohort_mut_count_tumour.pdf",
                   width = 6, height = 7)
  
  # plot mutation class distribution by cohort
  
  plots$P4 <-
    savePlotSafely(plot = plots$P4,
                   filename = "p4_cohort_mut_class.pdf",
                   width = 6, height = 5)
  
  # Case vs control error rates.
  
  plots$P5 <-
    savePlotSafely(plot = plots$P5,
                   filename = "p5_background_error_comparison.case_vs_control.pdf",
                   width = 6, height = 4)
  
  # Summary plots
  
  plots$P7 <-
    savePlotSafely(plot = plots$P7,
                   filename = "p7_background_error_rates.pdf",
                   width = 8, height = 4)
  
  plots$P9a <-
    savePlotSafely(plot = plots$P9a,
                   filename = "p9a_error_rate_comparison.both_strands.pdf",
                   width = 6, height = 3)
  
  plots$P9b <-
    savePlotSafely(plot = plots$P9b,
                   filename = "p9b_error_rate_comparison.locus_noise.pdf",
                   width = 6, height = 3)
  
  plots$P9c <-
    savePlotSafely(plot = plots$P9c,
                   filename = "p9c_error_rate_comparison.locus_noise_both_strands.pdf",
                   width = 6, height = 3)
  
  plots$P20 <-
    savePlotSafely(plot = plots$P20,
                   filename = "p20_filtering_comparison.pdf",
                   width = 6, height = 4)
  
  # Outlier suppression plots
  
  plots$P6 <-
    savePlotSafely(plot = plots$P6,
                   filename = "p6_os_data_retained.pdf",
                   width = 5, height = 6)
  
  plots$P8 <-
    savePlotSafely(plot = plots$P8,
                   filename = "p8_background_polishing.pdf",
                   width = 15, height = 12)
  
  # Tumour AF in observed and unobserved loci.
  
  plots$P10 <-
    savePlotSafely(plot = plots$P10,
                   filename = "p10_tumour_AF_for_observed_non_observed_loci.pdf",
                   width = 5, height = 4)
  
  # Fragment size per cohort, with different levels of error-suppression
  
  plots$P11 <-
    savePlotSafely(plot = plots$P11,
                   filename = "p11_size_comparison.pdf",
                   width = 6, height = 5)
  
  # Enrichment level. This can be NULL if there were no rows to create the plot from.
  
  if (is.null(plots$P12))
  {
    warning("There are no data points from which to create the enrichment level plot.")
  }
  else
  {
    plots$P12 <-
      savePlotSafely(plot = plots$P12,
                     filename = "p12_enrichment_ratios.pdf",
                     width = 6, height = 5)
  }
  
  # Receiver Operating Characteristic Plots
  
  plots$P13a <-
    savePlotSafely(plot = plots$P13a,
                   filename = "p13a_receiver_operating_characteristic.pdf",
                   width = 4, height = 3)
  
  plots$P13b <-
    savePlotSafely(plot = plots$P13b,
                   filename = "p13b_receiver_operating_characteristic.no_size.pdf",
                   width = 4, height = 3)
  
  # IR (depth) to IMAF plot
  
  if (is.null(plots$P14))
  {
    warning("There is no IF patient data from which to create the IF to IMAF plot.")
  }
  else
  {
    plots$P14 <-
      savePlotSafely(plot = plots$P14,
                     filename = "p14_IR_vs_IMAF.pdf",
                     width = 6, height = 4)
  }
  
  # Waterfall plot with detectable vs non detectable using dPCR
  
  if (is.null(plots$P15))
  {
    warning("There is no patient specific GLRT data from which to create the detectable/non detectable plot.")
  }
  else
  {
    plots$P15 <-
      savePlotSafely(plot = plots$P15,
                     filename = "p15_waterfall_IMAF.pdf",
                     width = 10, height = 5)
  }
  
  # Waterfall plot of cancer genomes
  
  if (is.null(plots$P16))
  {
    warning("There is no patient specific GLRT data from which to create the cancer genomes waterfall plot.")
  }
  else
  {
    plots$P16 <-
      savePlotSafely(plot = plots$P16,
                     filename = "p16_waterfall_cancer_genomes.pdf",
                     width = 10, height = 5)
  }
  
  # dPCR comparison plot
  
  if (is.null(plots$P17))
  {
    warning("There is no patient specific GLRT data from which to create the dPCR comparison plot.")
  }
  else
  {
    plots$P17 <-
      savePlotSafely(plot = plots$P17,
                     filename = "p17_dPCR_comparison.pdf",
                     width = 8, height = 5)
  }
  
  
  ## Render the INVAR analysis report.
  
  dir.create(str_c(getwd(), '/knitting'), showWarnings = FALSE)
  rmarkdown::render(input = str_c(Sys.getenv('INVAR_HOME'), '/R/5_analysis/analysis.Rmd'),
                    knit_root_dir = str_c(getwd(), '/knitting'),
                    intermediates_dir = str_c(getwd(), '/intermediates'),
                    output_format = rmarkdown::html_document(),
                    output_dir = getwd(),
                    output_file = str_c(scriptArgs$STUDY, "_invar2_analysis.html"))
}

# Launch it.

invisible(main(parseOptions()))
