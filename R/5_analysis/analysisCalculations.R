##
# Calculation functions for analysis.
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
    assert_that(!is.null(errorRatesTable), msg = "errorRatesTable cannot be null")
    assert_that(!is.null(layoutTable), msg = "layoutTable cannot be null")

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

annotatePatientSpecificGLRT <- function(patientSpecificGLRT, layoutTable, patientSummaryTable)
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
        mutate_at(vars(LS_FILTER, LOLLIPOP), as.factor) %>%
        left_join(patientSummaryTable, by = 'PATIENT') %>%
        mutate(CANCER_GENOMES_FRACTION = ifelse(DETECTED.WITH_SIZE, MUTATION_SUM / MUTATIONS, 2))

    patientSpecificGLRT.annotated
}

##
# From Emma's track_n_mut_v2.R script, translated for the new tibble structures.
#

mutationTracking <- function(mutationsTable, layoutTable, tumourMutationsTable, invarScoresTable)
{
    assert_that("PATIENT_SPECIFIC" %in% colnames(invarScoresTable),
                msg = "invarScoresTable lacking derived PATIENT_SPECIFIC column.")

    layoutTableTimepoint <- layoutTable %>%
        select(POOL, BARCODE, TIMEPOINT)

    tumourMutationTableSummary <- tumourMutationsTable %>%
        group_by(PATIENT) %>%
        summarise(INITIAL_MUTATIONS = n(), .groups = "drop")

    mutationsTable <- mutationsTable %>%
        mutate(UNIQUE_PATIENT_POS = str_c(UNIQUE_POS, UNIQUE_ALT, sep = '_'))

    # Only interested in patient specific rows from INVAR scores.
    # LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS are all true for the current
    # set up, so they don't need to be included in this table.

    invarScoresTable <- invarScoresTable %>%
        adjustInvarScores(layoutTable) %>%
        filter(PATIENT_SPECIFIC) %>%
        select(POOL, BARCODE, PATIENT, USING_SIZE, OUTLIER.PASS, DETECTION)

    mutationTracking <- mutationsTable %>%
        mutate(UNIQUE_IF_MUTANT_SPECIFIC = ifelse(MUTANT & PATIENT_SPECIFIC, UNIQUE_PATIENT_POS, NA),
               UNIQUE_IF_MUTANT_NON_SPECIFIC = ifelse(MUTANT & !PATIENT_SPECIFIC, UNIQUE_PATIENT_POS, NA)) %>%
        group_by(POOL, BARCODE, PATIENT, CASE_OR_CONTROL) %>%
        summarise(LOCUS_NOISE.PASS = sum(LOCUS_NOISE.PASS & PATIENT_SPECIFIC),
                  MUTANTS_PATIENT_SPECIFIC = n_distinct(UNIQUE_IF_MUTANT_SPECIFIC, na.rm = TRUE),
                  MUTANTS_NON_SPECIFIC = n_distinct(UNIQUE_IF_MUTANT_NON_SPECIFIC, na.rm = TRUE),
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
               MUTANTS_PATIENT_SPECIFIC, MUTANTS_NON_SPECIFIC,
               SPECIFIC_OUTLIER.PASS, SPECIFIC.PASS, NON_SPECIFIC.PASS,
               WITH_SIZE.OUTLIER_PASS, NO_SIZE.OUTLIER_PASS,
               WITH_SIZE.OUTLIER_FAIL, NO_SIZE.OUTLIER_FAIL) %>%
        arrange(POOL, BARCODE, PATIENT, TIMEPOINT, CASE_OR_CONTROL)

    mutationTracking
}
