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
        MUTATIONS_TABLE_FILE = str_c(base, 'mutation_table.with_outliers.rds'),
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
# Plot functions.
#

cohortMutationContextPlot <- function(contextMutationsTable, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    plot <- contextMutationsTable %>%
        ggplot(aes(x = as.character(TRINUCLEOTIDE), fill = MUTATION_CLASS)) +
            geom_bar() +
            theme_classic() +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14,face="bold"),
                  axis.text.x = element_text(angle = 90),
                  legend.position = "bottom") +
            labs(x = "Mutation context",
                 y = "Count",
                 title = str_c("Mutation context in ", study, " cohort"),
                 subtitle = "Unique patient specific mutations in polished filtered data")

    plot
}

cohortMutationAFByClassPlot <- function(contextMutationsTable, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    contextMutationsTable <- contextMutationsTable %>%
        filter(TUMOUR_AF > 0)

    plot <- contextMutationsTable %>%
        ggplot(aes(x = TUMOUR_AF, fill = MUTATION_CLASS)) +
            geom_histogram(binwidth = 0.02, position = "stack") +
            theme_classic() +
            labs(x = "AF tumour",
                 y = "Mutation count",
                 title = str_c("Mutation count in ", study, " tumour data"))+
            scale_fill_discrete(name = "Mutation")

    plot
}

mutationsPerPatientPlot <- function(patientSummaryTable, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    plot <- patientSummaryTable %>%
        ggplot(aes(x = PATIENT, y = MUTATIONS)) +
            geom_bar(stat = "identity") +
            theme_classic()+
            labs(x = "Patient",
                 y = "Tumour mutations",
                 title = str_c("Tumour mutation count in ", study, " cohort"))

    plot
}

mutationClassByCohortPlot <- function(contextMutationsTable, study)
{
    contextMutationsClassSummary <- contextMutationsTable %>%
        mutate(TOTAL_MUTATIONS = n()) %>%
        group_by(MUTATION_CLASS) %>%
        summarise(MUTATIONS = n(),
                  TOTAL_MUTATIONS = unique(TOTAL_MUTATIONS),
                  .groups = "drop") %>%
        mutate(FRACTION = MUTATIONS / TOTAL_MUTATIONS)

    plot <- contextMutationsClassSummary %>%
        ggplot(aes(x = MUTATION_CLASS, y = FRACTION)) +
            geom_bar(stat = "identity") +
            theme_classic() +
            theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14, face = "bold"),
                  axis.text.x = element_text(angle = 90),
                  legend.position = c(0.75, 0.8)) +
            labs(x = "Mutation class",
                 y = "Fraction",
                 title = str_c("Mutation class in ", study, " panel"),
                 subtitle = "Unique patient specific mutations in polished filtered data")

    plot
}


##
# From functions.R, get.error_rate_polishing_comparison
#

backgroundErrorRatesPlot <- function(errorRatesTable, study, familySize)
{
    plot <- errorRatesTable %>%
        filter(CASE_OR_CONTROL == 'case' & ERROR_RATE_TYPE == 'locus_noise.both_reads') %>%
        mutate_at(vars('TRINUCLEOTIDE'), as.factor) %>%
        ggplot(aes(x = reorder(TRINUCLEOTIDE, BACKGROUND_AF),  y = BACKGROUND_AF, colour = MUTATION_CLASS)) +
            geom_point() +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
                  axis.text=element_text(size=10),
                  axis.title=element_text(size=14)) +
            scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4)) +
            labs(x = "Trinucleotide context",
                 y = "Background error rate",
                 title = str_c("Background error rates for ", study, " fam_", familySize),
                 subtitle = "Using cases only") +
            guides(fill = 'none') +
            annotation_logticks(sides = "l")

    plot
}

errorRatePolishingComparisonPlot <- function(errorRatesTable, errorPolishingSetting)
{
    oneReadErrorRates <- errorRatesTable %>%
        filter(ERROR_RATE_TYPE == "one_read") %>%
        rename(BACKGROUND_AF.ONE_READ = BACKGROUND_AF)

    otherErrorRates <- errorRatesTable %>%
        filter(ERROR_RATE_TYPE == errorPolishingSetting) %>%
        rename(BACKGROUND_AF.OTHER = BACKGROUND_AF) %>%
        select(TRINUCLEOTIDE, ALT, CASE_OR_CONTROL, BACKGROUND_AF.OTHER)

    errorRatesComparison <- oneReadErrorRates %>%
        inner_join(otherErrorRates, by = c('TRINUCLEOTIDE', 'ALT', 'CASE_OR_CONTROL')) %>%
        mutate(RATIO = BACKGROUND_AF.OTHER / BACKGROUND_AF.ONE_READ)

    plotTitle <- errorPolishingSetting %>%
        str_replace_all('_', ' ') %>%
        str_replace_all('\\.', ' & ') %>%
        str_to_sentence()

    plot <- errorRatesComparison %>%
        filter(CASE_OR_CONTROL == 'case') %>%
        ggplot(aes(x = BACKGROUND_AF.ONE_READ, y = BACKGROUND_AF.OTHER, color = MUTATION_CLASS))+
        geom_point() +
        scale_y_log10(limits = c(0.5e-7, 1e-2)) +
        scale_x_log10(limits = c(0.5e-7, 1e-2)) +
        geom_abline(slope = 1, linetype = "dashed") +
        labs(x = "Error rate pre-filter",
             y = "Error rate post-filter",
             color = "Mutation class",
             title = plotTitle) +
        theme_classic()+
        guides(fill = 'none') +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))

    plot
}

filteringComparisonPlot <- function(errorRatesTable, study, familySize)
{
    plot <- errorRatesTable %>%
        filter(CASE_OR_CONTROL == 'case') %>%
        ggplot(aes(x = MUTATION_CLASS, y = BACKGROUND_AF, fill = ERROR_RATE_TYPE)) +
            geom_boxplot() +
            scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) +
            theme_classic() +
            scale_fill_discrete(name = "Background filtering")+
            labs(x = "Mutation class",
                 y = "Background AF",
                 fill = "Error rate type",
                 title = str_c("Background error filtering steps\n", study, ' fam_', familySize))

    plot
}

backgroundErrorCaseControlPlot <- function(errorRatesINV042, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    errorRatesINV042 %>%
        ggplot(aes(x = MUTATION_CLASS, y = BACKGROUND_AF, colour = CASE_OR_CONTROL)) +
            geom_boxplot()+
            scale_y_log10()+
            stat_compare_means(method = "t.test",
                               method.args = list(alternative = "greater"),
                               aes(label = str_c("p = ", ..p.format..))) +
            theme_bw() +
            labs(x = "Mutation class",
                 y = "Background error rate",
                 colour = "Case/control",
                 title = "Comparison of error rate estimates in case vs. control",
                 subtitle = str_c(study, " cohort")) +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14),
                  panel.grid.major = element_line(colour = alpha("black", 0.1))) +
            stat_compare_means(method = "t.test",
                               method.args = list(alternative = "greater"),
                               aes(label = str_c("p = ", ..p.format..)))
}

summaryCohortPlot <- function(mutationsTable, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    cohortSummary <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS &
               AF > 0 & AF < 0.25 & MUTATION_SUM < 10) %>%
        group_by(PATIENT_SPECIFIC) %>%
        summarise(PROPORTION = sum(OUTLIER.PASS) / n(),
                  TOTAL_ROWS = n(),
                  .groups = "drop") %>%
        mutate(COHORT = as.factor(ifelse(PATIENT_SPECIFIC, "Patient Specific", "Controls")),
               STUDY = study)

    plot <- cohortSummary %>%
        ggplot(aes(x = COHORT, y = PROPORTION, fill = STUDY, label = round(PROPORTION, digits = 2))) +
            geom_bar(position = "dodge", stat = "identity") +
            geom_text(vjust = 1.2) +
            labs(title = str_c(study, ": Effect of outlier suppression"),
                 fill = "Study",
                 x = "Cohort",
                 y = "Proportion of signal retained after filter") +
            theme_classic() +
            theme(legend.position = c(0.7, 0.7)) +
            scale_y_continuous(limits = c(0,1))

    plot
}

backgroundPolishingPlot <- function(mutationsTable, study, errorSuppression, outlierSuppression)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    # SAMPLE_NAME in original now has the patient mutation belongs to as part of it.
    # That might need to be done here too.

    plot <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS &
               AF > 0 & AF < 0.25 & MUTATION_SUM < 10) %>%
        mutate(COMBINED_SAMPLE_NAME = str_c(SAMPLE_NAME, " (", PATIENT_MUTATION_BELONGS_TO, ")"),
               PASS = ifelse(OUTLIER.PASS, "Yes", "No"),
               STUDY = study,
               COHORT = as.factor(ifelse(PATIENT_SPECIFIC, "Patient Specific", "Controls"))) %>%
        ggplot(aes(x = COMBINED_SAMPLE_NAME, y = AF, colour = PASS)) +
            geom_point() +
            theme_classic() +
            scale_y_log10() +
            facet_wrap(STUDY ~ COHORT, scales = "free_x", nrow = 3) +
            theme(axis.text.x=element_text(angle = 90)) +
            labs(title = str_c("Effect of outlier suppression ", study),
                 subtitle = str_c(errorSuppression, ", threshold = ", outlierSuppression),
                 colour = "Pass",
                 x = "Sample",
                 y = "Mutant allele fraction")

    plot
}

tumourAFInLociPlot <- function(mutationsTable, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    # recalculate average to exclude high level samples

    samplesToKeep <- mutationsTable %>%
        group_by(PATIENT, POOL_BARCODE) %>%
        summarise(AVERAGE = weighted.mean(AF, DP), .groups = "drop") %>%
        filter(AVERAGE <= 0.01) %>%
        distinct(POOL_BARCODE)

    mutationsTable.filtered <- mutationsTable %>%
        filter(OUTLIER.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS & LOCUS_NOISE.PASS &
               POOL_BARCODE %in% samplesToKeep$POOL_BARCODE) %>%
        mutate(IN_PLASMA = ifelse(AF > 0, "Yes", "No"),
               COHORT = as.factor(ifelse(PATIENT_SPECIFIC, "Patient Specific", "Non-Patient Specific")))

    plot <- mutationsTable.filtered %>%
        ggplot(aes(x = COHORT, y = TUMOUR_AF, fill = IN_PLASMA)) +
            geom_boxplot(outlier.colour = NA, notch = T) +
            theme_classic()+
            labs(x = "Data",
                 y = "Tumour allele fraction",
                 fill = "Observed in plasma",
                 title = str_c(study, ": Tumour AF of detected vs. non detected loci"),
                 subtitle = "Split by patient specific and controls") +
            theme(axis.text=element_text(size = 12),
                  axis.title=element_text(size = 14, face = "bold"),
                  panel.grid.major = element_line(colour = alpha("black", 0.1))) +
            stat_compare_means(method = "t.test",
                               method.args = list(alternative = "greater"),
                               aes(label = str_c("p = ", ..p.format..)))

    plot
}

## Fragment size per cohort, with different levels of error-suppression

fragmentSizePerCohortPlot <- function(sizeCharacterisationSummary, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    plot <- sizeCharacterisationSummary %>%
        filter(PATIENT_SPECIFIC) %>%
        mutate(STUDY = study,
               MUTANT_LABEL = as.factor(ifelse(MUTANT, "Yes", "No"))) %>%
        ggplot(aes(x = SIZE.ROUNDED, y = PROPORTION, fill = MUTANT_LABEL)) +
            geom_bar(stat = "identity", position = "dodge") +
            labs(x = "Fragment size in bp",
                 y = "Proportion",
                 fill = "Mutant",
                 title = str_c(study, ": Fragment size vs. error suppression")) +
            theme_classic() +
            theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14, face = "bold"),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 12, face = "bold")) +
            scale_x_continuous(limits = c(0, 500)) +
            facet_grid(STUDY~., scales = "free_y") +
            geom_vline(xintercept = c(166, 166*2), linetype = "dashed", alpha = 0.5) +
            annotate("text", x = 200, y = 0.15, label = "166bp", alpha = 0.5) +
            annotate("text", x = 370, y = 0.15, label = "332bp", alpha = 0.5)

    plot
}

enrichmentLevelPlot <- function(sizeCharacterisationSummary, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    joinColumns <- c("PATIENT_SPECIFIC", "CASE_OR_CONTROL", "SAMPLE_TYPE", "SIZE.ROUNDED")
    selectColumns <- c(joinColumns, "COUNT", "TOTAL", "PROPORTION")

    size.mutant = sizeCharacterisationSummary %>%
        filter(MUTANT) %>%
        select(all_of(selectColumns))

    size.wildType = sizeCharacterisationSummary %>%
        filter(!MUTANT) %>%
        select(all_of(selectColumns))

    enrichmentRatio <- size.wildType %>%
        left_join(size.mutant, by = joinColumns, suffix = c(".WT", ".M")) %>%
        group_by(PATIENT_SPECIFIC, CASE_OR_CONTROL) %>%
        mutate(TOTAL = sum(ifelse(is.na(PROPORTION.M), 0, COUNT.WT))) %>%
        ungroup() %>%
        mutate(PROPORTION.WT = COUNT.WT / TOTAL,
               RATIO = PROPORTION.M / PROPORTION.WT,
               RATIO.LOG2 = log2(RATIO))

    plot <- enrichmentRatio %>%
        mutate(STUDY = study) %>%
        filter(PATIENT_SPECIFIC & CASE_OR_CONTROL == 'case' & !is.na(RATIO.LOG2)) %>%
        ggplot(aes(x = SIZE.ROUNDED, y = RATIO.LOG2)) +
            geom_bar(stat = "identity", position = "dodge") +
            labs(x = "Fragment size in bp",
                 y = "Enrichment ratio - log2 ratio",
                 title = str_c(study, ": Enrichment ratio for ctDNA vs. error suppression")) +
            theme_classic() +
            theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14, face = "bold")) +
            scale_x_continuous(limits = c(0, 500)) +
            facet_grid(STUDY~.) +
            geom_vline(xintercept = c(166, 166*2), linetype = "dashed", alpha = 0.5) +
            annotate("text", x = 210, y = 5, label = "166bp", alpha = 0.5) +
            annotate("text", x = 380, y = 5, label = "332bp", alpha = 0.5)

    plot
}

receiverOperatingCharacteristicPlot <- function(invarScoresTable, layoutTable, withSizes, study, familySize)
{
    assert_that(is.logical(withSizes), msg = "withSizes must be a logical")

    cutoffName <- ifelse(withSizes, 'CUT_OFF.WITH_SIZE', 'CUT_OFF.NO_SIZE')

    adjustedScoresTable <- adjustInvarScores(invarScoresTable, layoutTable) %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & OUTLIER.PASS)

    scaledInvarResultsList <- adjustedScoresTable %>%
        scaleInvarScores()

    patientControlCutOffInfo <- scaledInvarResultsList[[cutoffName]]

    rocPatientControl <- patientControlCutOffInfo$ROC %>%
        mutate_at(vars(PATIENT_SPECIFIC), as.integer)

    patientControlSpecificity <- patientControlCutOffInfo$QUANTILE_SPECIFICITY
    patientControlCutOff <- patientControlCutOffInfo$INVAR_SCORE_THRESHOLD

    if (any(adjustedScoresTable$CASE_OR_CONTROL == "control_negative"))
    {
        # Data has healthy controls. Continue with analysis.

        healthyControlResultsList <- adjustedScoresTable %>%
            filter(PATIENT_SPECIFIC | CASE_OR_CONTROL == "control_negative") %>%
            mutate(CASE_OR_CONTROL = 'case') %>%
            scaleInvarScores()

        healthyCutOffInfo <- healthyControlResultsList[[cutoffName]]

        rocHealthyControl <- healthyCutOffInfo$ROC %>%
            mutate_at(vars(PATIENT_SPECIFIC), as.integer)

        healthyControlSpecificity <- healthyCutOffInfo$QUANTILE_SPECIFICITY
        healthyControlCutOff <- healthyCutOffInfo$INVAR_SCORE_THRESHOLD

        plot <- rocPatientControl %>%
            ggplot(aes(d = PATIENT_SPECIFIC, m = ADJUSTED_INVAR_SCORE)) +
                geom_roc(data = rocHealthyControl, cutoffs.at = healthyControlCutOff,
                         labels = FALSE, pointsize = 0.75, color = "red") + # to get rid of numbers on the ROC curve
                geom_roc(cutoffs.at = patientControlCutOff, labels = FALSE, pointsize = 0.75) + # to get rid of numbers on the ROC curve
                theme_classic() +
                labs(title = str_c(study, ", fam_", familySize),
                     x = "False positive fraction",
                     y = "True positive fraction",
                    subtitle = str_c("Specificity patients = ", patientControlSpecificity, "%\nSpecificity healthy = ", healthyControlSpecificity, "%"))
    }
    else
    {
        # This data does not have healthy controls, will compute the ROC curve from case data only

        plot <- rocPatientControl %>%
            ggplot(aes(d = PATIENT_SPECIFIC, m = ADJUSTED_INVAR_SCORE)) +
                geom_roc(cutoffs.at = patientControlCutOff, labels = FALSE, pointsize = 0.75) + # to get rid of numbers on the ROC curve
                theme_classic() +
                labs(title = str_c(study, ", fam_", familySize),
                     x = "False positive fraction",
                     y = "True positive fraction",
                     subtitle = str_c("Specificity patients = ", patientControlSpecificity, "%"))

    }

    plot
}

depthToIMAFPlot <- function(ifPatientData)
{
    tmp_20k <- tibble(X = c(5e3,5e3,20000,20000), Y = c(0, 2e-04, 5e-05, 0))

    tmp_66k <- tibble(X = c(20000,20000, 66666, 66666), Y = c(0,5e-05, 1.5e-5, 0))

    tmp_1m <- tibble(X = c(1e6,1e6, 5e6, 5e6), Y = c(1e-06, 3e-1, 3e-1, 2e-07))

    senseLine <- tibble(DP = c(5e3, 5e4, 5e5, 5e6), IMAF = 1 / DP)

    plot <- ifPatientData %>%
        mutate(LOGABLE_INVAR = ifelse(IMAF == 0, 1e-7, ADJUSTED_IMAF)) %>%
        mutate(DETECTED_LABEL = as.factor(ifelse(DETECTED.WITH_SIZE, "Yes", "No"))) %>%
        ggplot(aes(x = DP, y = LOGABLE_INVAR)) +
            geom_polygon(data = tmp_20k, aes(x = X, y = Y), alpha = 0.2, fill = "darkslategray4") +
            geom_polygon(data = tmp_66k, aes(x = X, y = Y), alpha = 0.2, fill = "darkslategray3") +
            geom_polygon(data = tmp_1m, aes(x = X, y = Y), alpha = 0.3, fill = "#F59C20") +
            geom_vline(xintercept = 20000, linetype = "dotted") +
            geom_vline(xintercept = 66666, linetype = "dotted") +
            geom_vline(xintercept = 1e6, linetype = "dotted") +
            geom_point(aes(colour = DETECTED_LABEL)) +
            geom_line(data = senseLine, aes(x = DP, y = IMAF), linetype = "longdash") +
            scale_x_log10(limits = c(5e3, 5e6)) +
            scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
                          labels = c("ND", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")) +
            theme_classic() +
            labs(x = "IR",
                 y = "IMAF",
                 colour = "Detected with size",
                 title = "PBCP IR vs IMAF plot")

    plot
}

detectableWaterfallPlot <- function(patientSpecificGLRT, layoutTable, study)
{
    layoutTable <- layoutTable %>%
        select(POOL, BARCODE, STUDY, INPUT_INTO_LIBRARY_NG, QUANTIFICATION_METHOD) %>%
        mutate_at(vars(INPUT_INTO_LIBRARY_NG), as.double)

    patientSpecificGLRT.annotated <- patientSpecificGLRT %>%
        left_join(layoutTable, by = c('POOL', 'BARCODE')) %>%
        mutate(NOT_DETECTABLE_DPCR = ADJUSTED_IMAF < 3.3 / INPUT_INTO_LIBRARY_NG,
               CTDNA_PLOTTING = ifelse(!DETECTED.WITH_SIZE, 1e-7, ifelse(DP < 20000, 1e-8, ADJUSTED_IMAF)),
               LS_FILTER = ifelse(DETECTED.WITH_SIZE | DP >= 20000 , "Pass", "Fail"),
               LOLLIPOP = ifelse(DETECTED.WITH_SIZE & NOT_DETECTABLE_DPCR, "non_dPCR", LS_FILTER)) %>%
        mutate_at(vars(LS_FILTER, LOLLIPOP), as.factor)

    plot <- patientSpecificGLRT.annotated %>%
        mutate(POOL_BARCODE = str_c(POOL, BARCODE, sep = " ")) %>%
        ggplot(aes(x = reorder(POOL_BARCODE, CTDNA_PLOTTING), y = (7.2 + log10(CTDNA_PLOTTING)))) +
            geom_bar(stat =  "identity", width = 0.5, colour = "white") +
            scale_fill_manual(values = "chocolate1") +
            geom_point(aes(shape = LOLLIPOP), size = 3) +
            scale_shape_manual(name = "",
                               limits = c("Fail", "non_dPCR"),
                               labels = c("Low sensitivity", "Non dPCR"),
                               values = c(1, 19)) +
            scale_y_continuous(breaks = c(-0.8, 0.2, 1.2, 2.2, 3.2, 4.2, 5.2, 6.2),
                               labels = c("ND", "ND", 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
            theme_classic() +
            labs(x = "Sample",
                 y = "IMAF",
                 title = str_c(study, ": All IMAFs")) +
            theme(axis.text.x=element_blank())

    plot
}

##
# Calculation functions.
#

##
# Complement the MUTATION_CLASS columns for reference alleles
# 'A' and 'G'. 'T' and 'C' remain unchanged.
#
convertComplementaryMutations <- function(backgroundErrorTable)
{
    complement <- function(sequence)
    {
        chartr('ATCG', 'TAGC', sequence)
    }

    complementary <- function(base)
    {
        base == 'A' | base == 'G'
    }

    forward <- backgroundErrorTable %>%
        filter(!complementary(REF))

    reverse <- backgroundErrorTable %>%
        filter(complementary(REF))

    reverse <- reverse %>%
        mutate(MUTATION_CLASS = complement(MUTATION_CLASS))

    bind_rows(forward, reverse) %>%
        arrange(TRINUCLEOTIDE, CASE_OR_CONTROL, REF, ALT)
}

calculateErrorRatesINV042 <- function(errorRatesTable, layoutTable)
{
    layoutTable <- layoutTable %>%
        select(POOL, BARCODE, CASE_OR_CONTROL)

    errorRatesTable <- errorRatesTable %>%
        left_join(layoutTable, by = c('POOL', 'BARCODE'))

    controls <- errorRatesTable %>%
        filter(str_detect(CASE_OR_CONTROL, "control"))

    cases <- errorRatesTable %>%
        filter(CASE_OR_CONTROL == "case") %>%
        slice_sample(n = nrow(controls), replace = TRUE)

    backgroundErrorTable <-
        bind_rows(cases, controls) %>%
        group_by(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL) %>%
        summarise(MUTATION_SUM = sum(MUTATION_SUM),
                  DP_SUM = sum(DP_SUM),
                  .groups = "drop") %>%
        group_by(TRINUCLEOTIDE, CASE_OR_CONTROL) %>%
        mutate(DP = sum(DP_SUM)) %>%
        ungroup() %>%
        filter(ALT != '.') %>%
        group_by(TRINUCLEOTIDE, CASE_OR_CONTROL, REF, ALT, DP) %>%
        summarise(MUTATION_SUM = sum(MUTATION_SUM),
                  .groups = "drop") %>%
        mutate(BACKGROUND_AF = MUTATION_SUM / DP) %>%
        mutate(ERROR_RATE_TYPE = 'one_strand', .after = 'CASE_OR_CONTROL') %>%
        mutate(MUTATION_CLASS = str_c(REF, ALT, sep = "/"))

    backgroundErrorTable %>%
        convertComplementaryMutations()
}

calculateErrorRateSummary <- function(errorRatesINV042)
{
    errorRatesINV042 %>%
        group_by(MUTATION_CLASS) %>%
        mutate(CLASS_ERROR_RATE = weighted.mean(BACKGROUND_AF)) %>%
        ungroup() %>%
        group_by(MUTATION_CLASS, TRINUCLEOTIDE) %>%
        mutate(CONTEXT_ERROR_RATE = weighted.mean(BACKGROUND_AF)) %>%
        ungroup() %>%
        summarise(MINIMUM_BY_CLASS = min(CLASS_ERROR_RATE),
                  MAXIMUM_BY_CLASS = max(CLASS_ERROR_RATE),
                  MINIMUM_BY_3BP = min(CONTEXT_ERROR_RATE),
                  MAXIMUM_BY_3BP = max(CONTEXT_ERROR_RATE))
}

calculateSizeCharacterisationSummary <- function(sizeCharacterisationTable, layoutTable, study, roundTo)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")
    assert_that(is.numeric(roundTo), msg = "roundTo must be a number.")

    studyInfo <- layoutTable %>%
        filter(STUDY == study) %>%
        distinct(STUDY, SAMPLE_TYPE)

    assert_that(nrow(studyInfo) == 1, msg = str_c("Different samples types in ", study))

    sizeCharacterisationTable %>%
        mutate(SIZE.ROUNDED = plyr::round_any(SIZE, accuracy = roundTo)) %>%
        group_by(PATIENT_SPECIFIC, CASE_OR_CONTROL, MUTANT, SIZE.ROUNDED) %>%
        summarise(COUNT = sum(COUNT), .groups = "drop_last") %>%
        mutate(TOTAL = sum(COUNT)) %>%
        ungroup() %>%
        mutate(PROPORTION = COUNT / TOTAL,
               SAMPLE_TYPE = studyInfo$SAMPLE_TYPE)
}

##
# Originally get.INVAR_score in TAPAS_functions.R
#

adjustInvarScores <- function(invarScoresTable, layoutTable)
{
    adjustThreshold <- function(row, conditions, adjustedScoresTable, scoreSpecificity)
    {
        condition <- slice(conditions, n = row)

        scores <- adjustedScoresTable %>%
            filter(USING_SIZE == condition$USING_SIZE &
                   LOCUS_NOISE.PASS == condition$LOCUS_NOISE.PASS &
                   BOTH_STRANDS.PASS == condition$BOTH_STRANDS.PASS &
                   OUTLIER.PASS == condition$OUTLIER.PASS) %>%
            filter(PATIENT_SPECIFIC | CONTAMINATION_RISK.PASS)

        scores.general <- scores %>%
            filter(!PATIENT_SPECIFIC) %>%
            arrange(ADJUSTED_INVAR_SCORE)

        percentilePosition = floor(nrow(scores.general) * scoreSpecificity)

        threshold <- scores.general$ADJUSTED_INVAR_SCORE[percentilePosition]

        scores %>%
            mutate(DETECTION = ADJUSTED_INVAR_SCORE > threshold)
    }

    layoutTable <- layoutTable %>%
        select(POOL, BARCODE, CASE_OR_CONTROL, SAMPLE_TYPE, DATA_TYPE, TIMEPOINT)

    assert_that(!any(invarScoresTable$DP < 0), msg = "Have negative DP.")

    # setting ctDNA level to 1/# molecules where p_mle (AF_P) < 1/# molecules (last mutate)

    adjustedScoresTable <- invarScoresTable %>%
        left_join(layoutTable, by = c('POOL', 'BARCODE')) %>%
        mutate(ADJUSTED_INVAR_SCORE = ifelse(MUTANT_READS_PRESENT, INVAR_SCORE, 0),
               ADJUSTED_IMAF = ifelse(MUTANT_READS_PRESENT, IMAF, 0)) %>%
        mutate(ADJUSTED_IMAF = ifelse(ADJUSTED_IMAF < 1 / DP, 1 / DP, ADJUSTED_IMAF))

    uniqueConditions <- invarScoresTable %>%
        distinct(USING_SIZE, LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS)

    adjustedList <- lapply(1:nrow(uniqueConditions), adjustThreshold,
                           uniqueConditions, adjustedScoresTable, scoreSpecificity = 0.95)

    adjustedScoresTable <-
        bind_rows(adjustedList) %>%
        arrange(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO,
                ITERATION, USING_SIZE, LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS)
}

##
# Originally scale_INVAR_size in functions.R
#

scaleInvarScores <- function(adjustedScoresTable, lowSensitivityThreshold = 20000)
{
    assert_that(is.numeric(lowSensitivityThreshold), msg = "lowSensitivityThreshold my be a number")

    specific <- adjustedScoresTable %>%
        filter(PATIENT_SPECIFIC)

    contaminationRiskSamples <- specific %>%
        filter(!USING_SIZE & ADJUSTED_IMAF > 0.01)

    nonSpecific <- adjustedScoresTable %>%
        filter(!PATIENT_SPECIFIC & CASE_OR_CONTROL == 'case' &
               !POOL_BARCODE %in% contaminationRiskSamples$POOL_BARCODE)

    joinColumns <- c('POOL', 'BARCODE', 'PATIENT', 'PATIENT_MUTATION_BELONGS_TO')
    joinSuffix <- c('.NO_SIZE', '.WITH_SIZE')

    nonSpecific.noSize <- nonSpecific %>%
        filter(!USING_SIZE) %>%
        select(all_of(joinColumns), ADJUSTED_INVAR_SCORE, DP)

    nonSpecific.withSize <- nonSpecific %>%
        filter(USING_SIZE) %>%
        select(all_of(joinColumns), ADJUSTED_INVAR_SCORE, IMAF, ADJUSTED_IMAF)

    beforeAfter.nonSpecific <- nonSpecific.noSize %>%
        left_join(nonSpecific.withSize, joinColumns, suffix = joinSuffix) %>%
        mutate(ADJUSTED_INVAR_SCORE.WITH_SIZE = ifelse(is.na(ADJUSTED_INVAR_SCORE.WITH_SIZE), 0, ADJUSTED_INVAR_SCORE.WITH_SIZE)) %>%
        filter(DP > lowSensitivityThreshold)

    joinColumns <- c('POOL', 'BARCODE')

    specific.noSize <- specific %>%
        filter(!USING_SIZE) %>%
        select(all_of(joinColumns), PATIENT, PATIENT_MUTATION_BELONGS_TO, ADJUSTED_INVAR_SCORE, DP, MUTATION_SUM, TIMEPOINT)

    specific.withSize <- specific %>%
        filter(USING_SIZE) %>%
        select(all_of(joinColumns), IMAF, ADJUSTED_INVAR_SCORE, ADJUSTED_IMAF)

    beforeAfter.specific <- specific.noSize %>%
        left_join(specific.withSize, joinColumns, suffix = joinSuffix)

    cutPointInfo.withSize <- cutPointGLRT(beforeAfter.specific, beforeAfter.nonSpecific, TRUE, lowSensitivityThreshold)

    cutPointInfo.noSize <- cutPointGLRT(beforeAfter.specific, beforeAfter.nonSpecific, FALSE, lowSensitivityThreshold)

    beforeAfter.nonSpecific <- beforeAfter.nonSpecific %>%
        mutate(DETECTED.NO_SIZE = ADJUSTED_INVAR_SCORE.NO_SIZE >= cutPointInfo.noSize$INVAR_SCORE_THRESHOLD,
               DETECTED.WITH_SIZE = ADJUSTED_INVAR_SCORE.WITH_SIZE >= cutPointInfo.withSize$INVAR_SCORE_THRESHOLD)

    beforeAfter.specific <- beforeAfter.specific %>%
        mutate(DETECTED.NO_SIZE = ADJUSTED_INVAR_SCORE.NO_SIZE >= cutPointInfo.noSize$INVAR_SCORE_THRESHOLD,
               DETECTED.WITH_SIZE = ADJUSTED_INVAR_SCORE.WITH_SIZE >= cutPointInfo.withSize$INVAR_SCORE_THRESHOLD)

    beforeAfter.nonSpecific.N <- nrow(beforeAfter.nonSpecific)

    beforeAfter.specific <- beforeAfter.specific %>%
        rowwise() %>%
        mutate(SPECIFICITY.WITH_SIZE =
                   sum(ADJUSTED_INVAR_SCORE.WITH_SIZE > beforeAfter.nonSpecific$ADJUSTED_INVAR_SCORE.WITH_SIZE) /
                   beforeAfter.nonSpecific.N) %>%
        mutate(SPECIFICITY.NO_SIZE =
                   sum(ADJUSTED_INVAR_SCORE.NO_SIZE > beforeAfter.nonSpecific$ADJUSTED_INVAR_SCORE.NO_SIZE) /
                   beforeAfter.nonSpecific.N) %>%
        ungroup()

    list(PATIENT_SPECIFIC = beforeAfter.specific,
         NON_SPECIFIC = beforeAfter.nonSpecific,
         CUT_OFF.WITH_SIZE = cutPointInfo.withSize,
         CUT_OFF.NO_SIZE = cutPointInfo.noSize)
}

##
# Originally cut_point_GLRT_no_size & cut_point_GLRT in functions.R
#

cutPointGLRT <- function(specificInvarScores, nonSpecificInvarScores, useSize, lowSensitivityThreshold)
{
    assert_that(is.logical(useSize), msg = "useSize must be a logical.")
    assert_that(is.numeric(lowSensitivityThreshold), msg = "lowSensitivityThreshold my be a number")

    invarScoreColumn <- ifelse(useSize, 'ADJUSTED_INVAR_SCORE.WITH_SIZE', 'ADJUSTED_INVAR_SCORE.NO_SIZE')
    useColumns <- c('POOL', 'BARCODE', 'PATIENT', 'PATIENT_MUTATION_BELONGS_TO', invarScoreColumn, 'DP')

    minimumSpecificDP <- min(specificInvarScores$DP)

    # low sensitivity threshold based on lowest patient sample

    roc <- bind_rows(specificInvarScores, nonSpecificInvarScores) %>%
        select(all_of(useColumns)) %>%
        mutate(PATIENT_SPECIFIC = PATIENT == PATIENT_MUTATION_BELONGS_TO) %>%
        rename(ADJUSTED_INVAR_SCORE = {{ invarScoreColumn }}) %>%
        filter(DP >= minimumSpecificDP)

    # optimal.cutpoints does not like a tibble!

    cutPointInfo <-
        OptimalCutpoints::optimal.cutpoints(data = as.data.frame(roc), X = 'ADJUSTED_INVAR_SCORE',
                                            status = "PATIENT_SPECIFIC", tag.healthy = FALSE,
                                            methods = "MaxSpSe", direction = "<")

    sizeCutOff <- cutPointInfo$MaxSpSe$Global$optimal.cutoff$cutoff
    maximumSpecificity <- max(cutPointInfo$MaxSpSe$Global$optimal.cutoff$Sp)

    scoresForQuantile <- nonSpecificInvarScores  %>%
        filter(DP >= minimumSpecificDP) %>%
        select({{ invarScoreColumn }})

    quantileResult <- quantile(scoresForQuantile[[1]], probs = maximumSpecificity)
    quantileLabel <- round(as.double(str_replace(names(quantileResult), '%', '')), digits = 1)
    invarScoreThreshold <- unname(quantileResult)

    # Original code changed the test to > 0 if the threshold is zero.
    # This minor adjustment saves the separate test and additional code by
    # setting the threshold to a tiny number above zero.

    list(SIZE_CUT_OFF = sizeCutOff,
         MAXIMUM_SPECIFICITY = maximumSpecificity,
         ROC = roc,
         QUANTILE_SPECIFICITY = quantileLabel,
         INVAR_SCORE_THRESHOLD = ifelse(invarScoreThreshold == 0, 1e-31, invarScoreThreshold))
}

getIFPatientData <- function(invarScoresTable, layoutTable, patientSummaryTable)
{
    adjustedScoresTable <- adjustInvarScores(invarScoresTable, layoutTable) %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & OUTLIER.PASS)

    scaledInvarResultsList <- adjustedScoresTable %>%
        scaleInvarScores()

    patientSpecificGLRT <- scaledInvarResultsList$PATIENT_SPECIFIC %>%
        arrange(POOL, BARCODE, PATIENT, PATIENT_MUTATION_BELONGS_TO)

    ifPatientData <- patientSpecificGLRT %>%
        left_join(patientSummaryTable, by = 'PATIENT') %>%
        mutate(UNIQUE_MOLECULES = DP / MUTATIONS,
               NG_ON_SEQ = UNIQUE_MOLECULES / 300,
               LOW_SENSITIVITY = DP < 20000 & !DETECTED.WITH_SIZE) %>%
        arrange(POOL, BARCODE, PATIENT, PATIENT_MUTATION_BELONGS_TO)

    thresholds <- as_tibble(c(0,20000,66666)) %>%
        rename(THRESHOLD = 1)

    thresholdEffects <- patientSpecificGLRT %>%
        crossing(thresholds) %>%
        group_by(THRESHOLD) %>%
        summarise(CASES = sum(DP >= THRESHOLD),
                  DETECTED = sum(DP >= THRESHOLD & DETECTED.WITH_SIZE),
                  NOT_DETECTED = sum(DP >= THRESHOLD & !DETECTED.WITH_SIZE),
                  DISCARDED = sum(DP < THRESHOLD),
                  .groups = "drop") %>%
        mutate(ALL_CASES = nrow(patientSpecificGLRT),
               ALL_DETECTED = sum(patientSpecificGLRT$DETECTED.WITH_SIZE),
               DETECTION_RATE_HARSH = DETECTED / CASES,
               DETECTION_RATE_SOFT = ALL_DETECTED / (ALL_DETECTED + NOT_DETECTED))

    list(PATIENT_SPECIFIC_GLRT = patientSpecificGLRT,
         IF_PATIENT_DATA = ifPatientData,
         THRESHOLD_EFFECTS = thresholdEffects)
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
        calculateErrorRatesINV042(offTargetErrorRatesList[['pre_filter']], layoutTable)

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

    ggsave(plot = errorRatePolishingComparisonPlot(errorRatesTable, "both_reads"),
           filename = "p9_error_rate_comparison.both_reads.pdf",
           width = 6, height = 3)

    ggsave(plot = errorRatePolishingComparisonPlot(errorRatesTable, "locus_noise"),
           filename = "p9_error_rate_comparison.locus_noise.pdf",
           width = 6, height = 3)

    ggsave(plot = errorRatePolishingComparisonPlot(errorRatesTable, "locus_noise.both_reads"),
           filename = "p9_error_rate_comparison.locus_noise_both_reads.pdf",
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
                                          study = scriptArgs$STUDY),
           filename = "p15_waterfall_IMAF.pdf",
           width = 10, height = 5)
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
