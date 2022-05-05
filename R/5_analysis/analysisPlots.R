##
# Plot functions for analysis.
#

# onTargetLocusErrorRatePlot isn't in analysis.R, but was from TAPAS_functions.R
# embedded in the annotate_with_locus_error_rate function. It was the case that
# this plot was a side effect of doing the calculations. It's better to have it
# here with the other plots and include it in the analysis report.

onTargetLocusErrorRatePlot <- function(errorRateTable, study, tapasSetting)
{
    nonZeroLoci <- errorRateTable %>%
        summarise(F = sum(BACKGROUND_AF > 0) / n())

    nonZeroLociPercentage <- signif(nonZeroLoci$F[1] * 100, digits = 3)

    locusNoiseFail <- errorRateTable %>%
        summarise(F = sum(!LOCUS_NOISE.PASS) / n())

    locusNoiseFailPercentage <- signif(locusNoiseFail$F[1] * 100, digits = 3)

    plot <- errorRateTable %>%
        ggplot(aes(x = BACKGROUND_AF, fill = LOCUS_NOISE.PASS)) +
        geom_histogram(bins = 100, position = "dodge") +
        scale_colour_discrete(name = "Locus noise pass") +
        scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1)) +
        theme_bw() +
        labs(x = "Background plasma AF across all samples",
             y = "Frequency",
             title = str_c(study, tapasSetting, "on-target, non-patient specific data", sep=" "),
             subtitle = str_c("Blacklisting of non-zero AF loci (representing ", nonZeroLociPercentage,
                              "% of data)\nBlacklisted loci (LOCUS_NOISE.FAIL) = ", locusNoiseFailPercentage,
                              "%\nsplit by COSMIC mutation status")) +
        facet_wrap(~COSMIC, scales = "free_y")

    plot
}


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
        ggplot(aes(x = reorder(PATIENT, MUTATIONS), y = MUTATIONS)) +
            geom_bar(stat = "identity") +
            theme_classic()+
            labs(x = "Patient",
                 y = "Tumour mutations",
                 title = str_c("Tumour mutation count in ", study, " cohort"))

    plot
}

inUsedMutationsPerPatientPlot <- function(patientSummaryTable, inputMutationsTable, study)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    # Making dataframe for plotting
    summaryTable <- left_join(patientSummaryTable, inputMutationsTable, by="PATIENT")
    df_in <- stack(summaryTable[,2:3])
    df_in=cbind(df_in ,patient=c(summaryTable[,1]))

    # Changing values for better understanding
    colnames(df_in) = c("N_MUTATIONS" , "MUT_TYPE"  ,   "PATIENT")
    # df_in[df_in=="MUTATIONS"]="MUTATIONS_POST_FILTERS"
    # df_in <- data.frame(lapply(df_in, function(x) {gsub("[!INPUT_MUTATIONS]", "MUTATIONS_POST_FILTERS", x)}))

    plot <- df_in %>%
        ggplot() +
        geom_bar(aes(x=PATIENT, y= N_MUTATIONS, fill=MUT_TYPE ), stat="identity",position = "dodge") +
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
        filter(CASE_OR_CONTROL == 'case' & ERROR_RATE_TYPE == 'locus_noise.both_strands') %>%
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
                 colour = "Mutation Class",
                 subtitle = "Using cases only") +
            guides(fill = 'none') +
            annotation_logticks(sides = "l")

    plot
}

errorRatePolishingComparisonPlot <- function(errorRatesTable, errorPolishingSetting)
{
    oneReadErrorRates <- errorRatesTable %>%
        filter(ERROR_RATE_TYPE == "prefilter") %>%
        rename(BACKGROUND_AF.PREFILTER = BACKGROUND_AF)

    otherErrorRates <- errorRatesTable %>%
        filter(ERROR_RATE_TYPE == errorPolishingSetting) %>%
        rename(BACKGROUND_AF.OTHER = BACKGROUND_AF) %>%
        select(TRINUCLEOTIDE, ALT, CASE_OR_CONTROL, BACKGROUND_AF.OTHER)

    errorRatesComparison <- oneReadErrorRates %>%
        inner_join(otherErrorRates, by = c('TRINUCLEOTIDE', 'ALT', 'CASE_OR_CONTROL')) %>%
        mutate(RATIO = BACKGROUND_AF.OTHER / BACKGROUND_AF.PREFILTER)

    plotTitle <- errorPolishingSetting %>%
        str_replace_all('_', ' ') %>%
        str_replace_all('\\.', ' & ') %>%
        str_to_sentence()

    plot <- errorRatesComparison %>%
        filter(CASE_OR_CONTROL == 'case') %>%
        ggplot(aes(x = BACKGROUND_AF.PREFILTER, y = BACKGROUND_AF.OTHER, color = MUTATION_CLASS))+
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
        mutate(ERROR_RATE_TYPE = factor(ERROR_RATE_TYPE, levels = c('prefilter', 'both_strands', 'locus_noise', 'locus_noise.both_strands'))) %>%
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

backgroundPolishingPlot <- function(mutationsTable, layoutTable, study, errorSuppression, outlierSuppression)
{
    assert_that(is.character(study), msg = "Study is expected to be a string")

    layoutTable <- layoutTable %>%
        select(SAMPLE_ID, SAMPLE_NAME)

    # SAMPLE_NAME in original now has the patient mutation belongs to as part of it.
    # That might need to be done here too.

    plot <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS &
               AF > 0 & AF < 0.25 & MUTATION_SUM < 10) %>%
        left_join(layoutTable, by = 'SAMPLE_ID') %>%
        mutate(COMBINED_SAMPLE_NAME = str_c(SAMPLE_NAME, " (", PATIENT_MUTATION_BELONGS_TO, ")"),
               PASS = ifelse(OUTLIER.PASS, "Yes", "No"),
               STUDY = study,
               COHORT = as.factor(ifelse(PATIENT_SPECIFIC, "Patient Specific", "Non Specific"))) %>%
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
        group_by(PATIENT, SAMPLE_ID) %>%
        summarise(AVERAGE = weighted.mean(AF, DP), .groups = "drop") %>%
        filter(AVERAGE <= 0.01) %>%
        distinct(SAMPLE_ID)

    mutationsTable.filtered <- mutationsTable %>%
        filter(OUTLIER.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS & LOCUS_NOISE.PASS &
               SAMPLE_ID %in% samplesToKeep$SAMPLE_ID) %>%
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

    enrichmentRatio.filtered <- enrichmentRatio %>%
        filter(PATIENT_SPECIFIC & CASE_OR_CONTROL == 'case' & !is.na(RATIO.LOG2))

    plot <- NULL
    if (nrow(enrichmentRatio.filtered) > 0)
    {
        plot <- enrichmentRatio.filtered %>%
            mutate(STUDY = study) %>%
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
    }

    plot
}

receiverOperatingCharacteristicPlot <- function(invarScoresTable, layoutTable, withSizes, study, familySize, scoreSpecificity)
{
    assert_that(is.logical(withSizes), msg = "withSizes must be a logical")

    cutoffName <- ifelse(withSizes, 'CUT_OFF.WITH_SIZE', 'CUT_OFF.NO_SIZE')

    adjustedScoresTable <- adjustInvarScores(invarScoresTable, layoutTable, scoreSpecificity) %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & OUTLIER.PASS)

    scaledInvarResultsList <-
        tryCatch(
        {
            adjustedScoresTable %>%
                scaleInvarScores()
        },
        error = function(cond)
        {
            warning(geterrmessage())
            NULL
        })

    if (is.null(scaledInvarResultsList))
    {
        return(NULL)
    }

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
    assert_that(!is.null(ifPatientData), msg = "ifPatientData is null")

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

detectableWaterfallPlot <- function(annotatedPatientSpecificGLRT, study)
{
    plot <- annotatedPatientSpecificGLRT %>%
        ggplot(aes(x = reorder(SAMPLE_ID, CTDNA_PLOTTING), y = (7.2 + log10(CTDNA_PLOTTING)))) +
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

cancerGenomesWaterfallPlot <- function(annotatedPatientSpecificGLRT, study)
{
    plot <- annotatedPatientSpecificGLRT %>%
        ggplot(aes(x = reorder(SAMPLE_ID, CTDNA_PLOTTING), y = CANCER_GENOMES_FRACTION)) +
           geom_bar(stat =  "identity", width = 0.7, colour = "white") +
           scale_fill_manual(values = "chocolate1") +
           scale_y_log10(breaks = c(1e-4,1e-3,1e-2,1e-1,1,2,10,100),
                         labels = c(1e-4,1e-3,1e-2,1e-1,1,"ND",10,100)) +
           theme_classic() +
           labs(x = "Samples, ordered by IMAF",
                y = "Cancer genomes in sample") +
           theme(axis.text.x=element_blank())

    plot
}

dpcrComparisionPlot <- function(annotatedPatientSpecificGLRT, study)
{
    # 1ng = 300AC or 150GE
    # we are assuming 50% library efficiency
    # theoretical sensitivity with one locus assay would be 1/AC_input

    ## whats the proportion of detected samples that would not be detected with dPCR

    undetected <- annotatedPatientSpecificGLRT %>%
        filter(DETECTED.WITH_SIZE) %>%
        summarise(UNDETECTED = sum(NOT_DETECTABLE_DPCR) / n()) %>%
        as.double()

    cutoffLine <- tibble(NG_INPUT = 1:80) %>%
        mutate(CUTOFF = 3.3 / (NG_INPUT * 300))

    plot <- annotatedPatientSpecificGLRT %>%
        filter(ADJUSTED_IMAF > 0 & DETECTED.WITH_SIZE) %>%
        ggplot(aes(x = INPUT_INTO_LIBRARY_NG / 300, y = ADJUSTED_IMAF))+
        geom_line(data = cutoffLine, aes(x = NG_INPUT, y = CUTOFF))+
        geom_point(aes(color = NOT_DETECTABLE_DPCR)) +
        theme_classic() +
        labs(x = "Library input (ng)",
             y = "IMAF",
             title = str_c(study, ": ctDNA level of detected samples"),
             subtitle = str_c("line = 3.3/AC input; ", undetected * 100, "% of samples undetected")) +
        scale_y_log10(breaks=c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
        xlim(0,80) +
        scale_color_manual(breaks = c("FALSE", "TRUE"),
                           labels = c("Detectable", "Not detectable"),
                           name = "Detection with dPCR",
                           values = c("dodgerblue", "brown2"))

    plot
}


##
# Function to render plots and trap any errors that are raised while
# doing so, printing a warning instead. Returns the plot object passed
# in if there is no error, nor NULL if there was an error rendering
# the plot.
#

savePlotSafely <- function(plot, filename, width, height)
{
    if (is.null(plot))
    {
        return(NULL)
    }

    tryCatch(
    {
        ggsave(plot = plot, filename = filename,
               width = width, height = height)
        plot
    },
    error = function(cond)
    {
        warning("Could not render or save '", filename, "': ", geterrmessage())
        NULL
    })
}
