suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
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
        make_option(c("--family-size"), type="character", metavar="string",
                    dest="FAMILY_SIZE", help="The family size string",
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
        LAYOUT_FILE = str_c(testhome, 'source_files/combined.SLX_table_with_controls_031220.csv'),
        STUDY = 'PARADIGM',
        ERROR_SUPPRESSION = 'f0.9_s2',
        FAMILY_SIZE = 'fam_2',
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

cohortMutationAFByClassPlot <- function(contextMutationsTable, layoutTable, study)
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
        group_by(STUDY) %>%
        mutate(TOTAL_MUTATIONS = n()) %>%
        group_by(STUDY, MUTATION_CLASS) %>%
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
                 title = str_c("Background error rates for", study, familySize, sep = ' '),
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
                 title = str_c("Background error filtering steps\n", study, ' ', familySize))

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
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS & CONTAMINATION_RISK.PASS &
               AF > 0 & AF < 0.25 & MUTATION_SUM < 10) %>%
        group_by(STUDY, PATIENT_SPECIFIC) %>%
        summarise(PROPORTION = sum(OUTLIER.PASS) / n(),
                  TOTAL_ROWS = n(),
                  .groups = "drop") %>%
        mutate(COHORT = as.factor(ifelse(PATIENT_SPECIFIC, "Patient Specific", "Controls")))

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
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS & CONTAMINATION_RISK.PASS &
               AF > 0 & AF < 0.25 & MUTATION_SUM < 10) %>%
        mutate(COMBINED_SAMPLE_NAME = str_c(SAMPLE_NAME, " (", PATIENT_MUTATION_BELONGS_TO, ")"),
               PASS = ifelse(OUTLIER.PASS, "Yes", "No"),
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
        filter(OUTLIER.PASS & BOTH_STRANDS & CONTAMINATION_RISK.PASS & LOCUS_NOISE.PASS &
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

fragmentSizePerCohortPlot <- function(sizeCharacterisationTable, layoutTable, study)
{
    studyInfo <- layoutTable %>%
        filter(STUDY == study) %>%
        distinct(STUDY, SAMPLE_TYPE)

    assert_that(nrow(studyInfo) == 1, msg = str_c("Different samples types in ", study))

    roundTo <- 5

    summary <- sizeCharacterisationTable %>%
        mutate(STUDY = study,
               SIZE.ROUNDED = plyr::round_any(SIZE, accuracy = roundTo)) %>%
        group_by(STUDY, PATIENT_SPECIFIC, CASE_OR_CONTROL, MUTANT, SIZE.ROUNDED) %>%
        summarise(COUNT = sum(COUNT), .groups = "drop_last") %>%
        mutate(TOTAL = sum(COUNT)) %>%
        ungroup() %>%
        mutate(PROPORTION = COUNT / TOTAL,
               SAMPLE_TYPE = studyInfo$SAMPLE_TYPE)

    plot <- summary %>%
        filter(PATIENT_SPECIFIC) %>%
        mutate(MUTANT_LABEL = as.factor(ifelse(MUTANT, "Yes", "No"))) %>%
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

errorRateSummary <- function(errorRatesINV042)
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

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    layoutTable <-
        loadLayoutTable(scriptArgs$LAYOUT_FILE) %>%
        select(STUDY, POOL, BARCODE, SAMPLE_TYPE, CASE_OR_CONTROL, TIMEPOINT)

    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    sizeCharacterisationTable <- readRDS(scriptArgs$SIZE_CHARACTERISATION_FILE)

    invarScoresTable <- readRDS(scriptArgs$INVAR_SCORES_FILE)

    errorRatesTable <- readRDS(scriptArgs$ERROR_RATES_FILE) %>%
        mutate(MUTATION_CLASS = str_c(REF, ALT, sep = '/'))

    offTargetErrorRatesList <- readRDS(scriptArgs$OFF_TARGET_ERROR_RATES_FILE)

    # Manipulation and further calculations.

    contextMutationsTable <- mutationsTable %>%
        filter(PATIENT_SPECIFIC & LOCUS_NOISE.PASS & BOTH_STRANDS & OUTLIER.PASS) %>%
        group_by(PATIENT, UNIQUE_POS) %>%
        slice_head(n = 1) %>%
        ungroup()

    patientSummaryTable <- contextMutationsTable %>%
        group_by(STUDY, PATIENT) %>%
        summarise(MUTATIONS = n_distinct(STUDY, PATIENT, UNIQUE_POS), .groups = "drop") %>%
        arrange(PATIENT, MUTATIONS)

    patientSummaryTable %>%
        write_csv("tumour_mutation_per_patient.csv")

    errorRatesINV042 <-
        calculateErrorRatesINV042(offTargetErrorRatesList[['pre_filter']], layoutTable)

    errorRateSummary(errorRatesINV042) %>%
        write_csv("error_rate_summary.csv")

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


    ggsave(plot = fragmentSizePerCohortPlot(sizeCharacterisationTable, layoutTable, study = scriptArgs$STUDY),
           filename = "p11_size_comparison.pdf",
           width = 6, height = 5)
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
