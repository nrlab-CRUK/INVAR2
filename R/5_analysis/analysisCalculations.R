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
    assert_that(all(str_length(base) == 1), msg = "Have one or more strings in vector whose length is not 1")

    # See https://builtin.com/data-science/and-in-r for "&" vs "&&".
    base == 'A' | base == 'G'
  }

  forward <- backgroundErrorTable %>%
    filter(!complementary(REF))

  reverse <- backgroundErrorTable %>%
    filter(complementary(REF)) %>%
    mutate(MUTATION_CLASS = complement(MUTATION_CLASS))

  bind_rows(forward, reverse) %>%
    arrange(TRINUCLEOTIDE, CASE_OR_CONTROL, REF, ALT)
}

calculateErrorRatesINV042 <- function(errorRatesTable, layoutTable)
{
  assert_that(!is.null(errorRatesTable), msg = "errorRatesTable cannot be null")
  assert_that(!is.null(layoutTable), msg = "layoutTable cannot be null")

  layoutTable <- layoutTable %>%
    select(SAMPLE_ID, CASE_OR_CONTROL)

  errorRatesTable <- errorRatesTable %>%
    left_join(layoutTable, by = 'SAMPLE_ID')

  controls <- errorRatesTable %>%
    filter(str_detect(CASE_OR_CONTROL, "control"))

  cases <- errorRatesTable %>%
    filter(CASE_OR_CONTROL == "case")

  if (nrow(cases) > 0)
  {
    cases <- cases %>%
      slice_sample(n = nrow(controls), replace = TRUE)
  }

  backgroundErrorTable <-
    bind_rows(cases, controls) %>%
    group_by(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL) %>%
    summarise(MUTATED_READS_PER_LOCI = sum(MUTATED_READS_PER_LOCI),
              DP_SUM = sum(DP_SUM),
              .groups = "drop") %>%
    group_by(TRINUCLEOTIDE, CASE_OR_CONTROL) %>%
    mutate(DP = sum(DP_SUM)) %>%
    ungroup() %>%
    filter(ALT != '.') %>%
    group_by(TRINUCLEOTIDE, CASE_OR_CONTROL, REF, ALT, DP) %>%
    summarise(MUTATED_READS_PER_LOCI = sum(MUTATED_READS_PER_LOCI),
              .groups = "drop") %>%
    mutate(BACKGROUND_AF = MUTATED_READS_PER_LOCI / DP) %>%
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

adjustInvarScores <- function(invarScoresTable, layoutTable, scoreSpecificity)
{
  assert_that(is.numeric(scoreSpecificity), msg = "scoreSpecificity must be a number")

  adjustThreshold <- function(row, conditions, adjustedScoresTable, scoreSpecificity)
  {
    condition <- slice(conditions, row)

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
    select(SAMPLE_ID, CASE_OR_CONTROL, SAMPLE_TYPE, TIMEPOINT)

  assert_that(!any(invarScoresTable$DP < 0), msg = "Have negative DP.")

  # setting ctDNA level to 1/# molecules where p_mle (AF_P) < 1/# molecules (last mutate)

  adjustedScoresTable <- invarScoresTable %>%
    left_join(layoutTable, by = 'SAMPLE_ID') %>%
    mutate(ADJUSTED_INVAR_SCORE = ifelse(MUTANT_READS_PRESENT, INVAR_SCORE, 0),
           ADJUSTED_IMAF = ifelse(MUTANT_READS_PRESENT, IMAF, 0)) %>%
    mutate(ADJUSTED_IMAF = ifelse(ADJUSTED_IMAF < 1 / DP, 1 / DP, ADJUSTED_IMAF))

  uniqueConditions <- invarScoresTable %>%
    distinct(USING_SIZE, LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS)

  adjustedList <- lapply(1:nrow(uniqueConditions), adjustThreshold,
                         uniqueConditions, adjustedScoresTable, scoreSpecificity)

  adjustedScoresTable <-
    bind_rows(adjustedList) %>%
    arrange(SAMPLE_ID, PATIENT_MUTATION_BELONGS_TO,
            ITERATION, USING_SIZE, LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS)
}

##
# Originally scale_INVAR_size in functions.R
#

scaleInvarScores <- function(adjustedScoresTable, minInformativeReads, maxBackgroundAlleleFreq)
{
  assert_that(is.numeric(minInformativeReads), msg = "minInformativeReads must be a number")
  assert_that(is.numeric(maxBackgroundAlleleFreq), msg = "maxBackgroundAlleleFreq must be a number")

  specific <- adjustedScoresTable %>%
    filter(PATIENT_SPECIFIC)

  contaminationRiskSamples <- specific %>%
    filter(!USING_SIZE & ADJUSTED_IMAF > maxBackgroundAlleleFreq)

  nonSpecific <- adjustedScoresTable %>%
    filter(!PATIENT_SPECIFIC & CASE_OR_CONTROL == 'case' &
             !SAMPLE_ID %in% contaminationRiskSamples$SAMPLE_ID)

  joinColumns <- c('SAMPLE_ID', 'PATIENT', 'PATIENT_MUTATION_BELONGS_TO', 'ITERATION')
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
    filter(DP > minInformativeReads)

  joinColumns <- 'SAMPLE_ID'

  specific.noSize <- specific %>%
    filter(!USING_SIZE) %>%
    select(all_of(joinColumns), PATIENT, PATIENT_MUTATION_BELONGS_TO, ADJUSTED_INVAR_SCORE, DP, MUTATED_READS_PER_LOCI, TIMEPOINT)

  specific.withSize <- specific %>%
    filter(USING_SIZE) %>%
    select(all_of(joinColumns), IMAF, ADJUSTED_INVAR_SCORE, ADJUSTED_IMAF)

  beforeAfter.specific <- specific.noSize %>%
    left_join(specific.withSize, joinColumns, suffix = joinSuffix)

  cutPointInfo.withSize <- cutPointGLRT(beforeAfter.specific, beforeAfter.nonSpecific, TRUE)

  cutPointInfo.noSize <- cutPointGLRT(beforeAfter.specific, beforeAfter.nonSpecific, FALSE)

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

cutPointGLRT <- function(specificInvarScores, nonSpecificInvarScores, useSize)
{
  assert_that(is.logical(useSize), msg = "useSize must be a logical.")

  invarScoreColumn <- ifelse(useSize, 'ADJUSTED_INVAR_SCORE.WITH_SIZE', 'ADJUSTED_INVAR_SCORE.NO_SIZE')
  useColumns <- c('SAMPLE_ID', 'PATIENT', 'PATIENT_MUTATION_BELONGS_TO', invarScoreColumn, 'DP')

  minimumSpecificDP <- ifelse(nrow(specificInvarScores) == 0, 0, min(specificInvarScores$DP))

  # low sensitivity threshold based on lowest patient sample

  roc <- bind_rows(specificInvarScores, nonSpecificInvarScores) %>%
    select(all_of(useColumns)) %>%
    mutate(PATIENT_SPECIFIC = PATIENT == PATIENT_MUTATION_BELONGS_TO) %>%
    rename(ADJUSTED_INVAR_SCORE = {{ invarScoreColumn }}) %>%
    filter(DP >= minimumSpecificDP) # depth can't be lower then lowest patient specific depth

  # optimal.cutpoints does not like a tibble!

  cutPointInfo <-
    tryCatch(
      {
        OptimalCutpoints::optimal.cutpoints(data = as.data.frame(roc), X = 'ADJUSTED_INVAR_SCORE',
                                            status = "PATIENT_SPECIFIC", tag.healthy = FALSE,
                                            methods = "MaxSpSe", direction = "<")
      },
      error = function(cond)
      {
        stop("OptimalCutpoints::optimal.cutpoints failed: ", geterrmessage())
      })

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

getIFPatientData <- function(invarScoresTable, layoutTable, patientSummaryTable,
                             scoreSpecificity, minInformativeReads, maxBackgroundAlleleFreq)
{
  assert_that(is.numeric(scoreSpecificity), msg = "scoreSpecificity must be a number")
  assert_that(is.numeric(minInformativeReads), msg = "minInformativeReads must be a number")
  assert_that(is.numeric(maxBackgroundAlleleFreq), msg = "maxBackgroundAlleleFreq must be a number")

  adjustedScoresTable <- adjustInvarScores(invarScoresTable, layoutTable, scoreSpecificity) %>%
    filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & OUTLIER.PASS)

  scaledInvarResultsList <-
    tryCatch(
      {
        scaleInvarScores(adjustedScoresTable, minInformativeReads, maxBackgroundAlleleFreq)
      },
      error = function(cond)
      {
        stop("Patient IF data unavailable because of error: ", geterrmessage())
      })

  patientSpecificGLRT <- scaledInvarResultsList$PATIENT_SPECIFIC %>%
    arrange(SAMPLE_ID, PATIENT, PATIENT_MUTATION_BELONGS_TO)

  ifPatientData <- patientSpecificGLRT %>%
    left_join(patientSummaryTable, by = 'PATIENT') %>%
    mutate(UNIQUE_MOLECULES = DP / MUTATIONS,
           NG_ON_SEQ = UNIQUE_MOLECULES / 300,
           LOW_SENSITIVITY = DP < minInformativeReads & !DETECTED.WITH_SIZE) %>%
    arrange(SAMPLE_ID, PATIENT, PATIENT_MUTATION_BELONGS_TO)

  thresholds <- as_tibble(c(0, minInformativeReads, 66666)) %>%
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

annotatePatientSpecificGLRT <- function(patientSpecificGLRT, layoutTable, patientSummaryTable, minInformativeReads)
{
  assert_that(is.numeric(minInformativeReads), msg = "minInformativeReads must be a number")

  layoutTable <- layoutTable %>%
    select(SAMPLE_ID, STUDY, INPUT_INTO_LIBRARY_NG) %>%
    mutate_at(vars(INPUT_INTO_LIBRARY_NG), as.double)

  patientSpecificGLRT.annotated <- patientSpecificGLRT %>%
    left_join(layoutTable, by = 'SAMPLE_ID') %>%
    mutate(NOT_DETECTABLE_DPCR = ADJUSTED_IMAF < 3.3 / INPUT_INTO_LIBRARY_NG,
           CTDNA_PLOTTING = ifelse(!DETECTED.WITH_SIZE, 1e-7, ifelse(DP < minInformativeReads, 1e-8, ADJUSTED_IMAF)),
           LS_FILTER = ifelse(DETECTED.WITH_SIZE | DP >= minInformativeReads , "Pass", "Fail"),
           LOLLIPOP = ifelse(DETECTED.WITH_SIZE & NOT_DETECTABLE_DPCR, "non_dPCR", LS_FILTER)) %>%
    mutate_at(vars(LS_FILTER, LOLLIPOP), as.factor) %>%
    left_join(patientSummaryTable, by = 'PATIENT') %>%
    mutate(CANCER_GENOMES_FRACTION = ifelse(DETECTED.WITH_SIZE, MUTATED_READS_PER_LOCI / MUTATIONS, 2))

  patientSpecificGLRT.annotated
}

##
# From Emma's track_n_mut_v2.R script, translated for the new tibble structures.
#

mutationTracking <- function(mutationsTable, layoutTable, tumourMutationsTable, invarScoresTable, scoreSpecificity)
{
  assert_that("PATIENT_SPECIFIC" %in% colnames(invarScoresTable),
              msg = "invarScoresTable lacking derived PATIENT_SPECIFIC column.")

  layoutTableTimepoint <- layoutTable %>%
    select(SAMPLE_ID, TIMEPOINT)

  tumourMutationTableSummary <- tumourMutationsTable %>%
    group_by(PATIENT) %>%
    summarise(INITIAL_MUTATIONS = n(), .groups = "drop")

  mutationsTable <- mutationsTable %>%
    mutate(UNIQUE_PATIENT_POS = str_c(UNIQUE_POS, UNIQUE_ALT, sep = '_'),
           UNIQUE_IF_MUTANT_SPECIFIC = ifelse(MUTANT & PATIENT_SPECIFIC, UNIQUE_PATIENT_POS, NA),
           UNIQUE_IF_MUTANT_NON_SPECIFIC = ifelse(MUTANT & !PATIENT_SPECIFIC, UNIQUE_PATIENT_POS, NA),
           # UNIQUE_IF_MUTANT_CASE_OR_CONTROL = ifelse(CASE_OR_CONTROL == 'case', !is.na(UNIQUE_IF_MUTANT_SPECIFIC), !is.na(UNIQUE_IF_MUTANT_NON_SPECIFIC)),
           PASS_ALL = LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS & OUTLIER.PASS & MUTANT, # MUTATED_READS_PER_LOCI > 0
           ALL_IR = LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS & OUTLIER.PASS & PATIENT==PATIENT_MUTATION_BELONGS_TO,
           ALL_IR_CONTROL = LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS & OUTLIER.PASS)

  # Only interested in patient specific rows from INVAR scores.
  # LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & CONTAMINATION_RISK.PASS are all true for the current
  # set up, so they don't need to be included in this table.

  invarScoresTable <- invarScoresTable %>%
    adjustInvarScores(layoutTable, scoreSpecificity) %>%
    filter(PATIENT_SPECIFIC) %>%
    select(SAMPLE_ID, PATIENT, USING_SIZE, OUTLIER.PASS, DETECTION)

  mutationTracking <- mutationsTable %>%
    group_by(SAMPLE_ID, PATIENT, CASE_OR_CONTROL) %>%
    summarise(RAW_N_READS = sum(REF_F + REF_R + ALT_F + ALT_R),
              N_LOCI_MUTATED_PTSPEC = n_distinct(UNIQUE_IF_MUTANT_SPECIFIC, na.rm = TRUE),
              N_LOCI_MUTATED_NON_PTSPEC = n_distinct(UNIQUE_IF_MUTANT_NON_SPECIFIC, na.rm = TRUE),
              #
              N_READS_MUTATED_PTSPEC = sum(ALT_F[MUTANT & PATIENT_SPECIFIC] + ALT_R[MUTANT & PATIENT_SPECIFIC]), # Note this is number of reads NOT fragments (N_FRAGMENTS = 0.5 * N_READS)
              N_READS_MUTATED_NON_PTSPEC = sum(ALT_F[MUTANT & !PATIENT_SPECIFIC] + ALT_R[MUTANT & !PATIENT_SPECIFIC]),
              #
              N_LOCI_MUTATED_PTSPEC_LNP = n_distinct(UNIQUE_IF_MUTANT_SPECIFIC[LOCUS_NOISE.PASS], na.rm = TRUE),
              N_LOCI_MUTATED_NON_PTSPEC_LNP = n_distinct(UNIQUE_IF_MUTANT_NON_SPECIFIC[LOCUS_NOISE.PASS], na.rm = TRUE),
              # Adding the alt and reverse number of reads if UNIQUE_IF_MUTANT_SPECIFIC is not NA and LNP==TRUE
              N_READS_MUTATED_PTSPEC_LNP = sum(ALT_F[MUTANT & PATIENT_SPECIFIC & LOCUS_NOISE.PASS] + ALT_R[MUTANT & PATIENT_SPECIFIC & LOCUS_NOISE.PASS]), # To check if its the same as two lines above
              N_READS_MUTATED_NON_PTSPEC_LNP = sum(ALT_F[MUTANT & !PATIENT_SPECIFIC & LOCUS_NOISE.PASS] + ALT_R[MUTANT & !PATIENT_SPECIFIC & LOCUS_NOISE.PASS]),
              #
              N_LOCI_MUTATED_PTSPEC_ALL_FILTERS = n_distinct(UNIQUE_IF_MUTANT_SPECIFIC[PASS_ALL], na.rm = TRUE),
              N_LOCI_MUTATED_NON_PTSPEC_ALL_FILTERS = n_distinct(UNIQUE_IF_MUTANT_NON_SPECIFIC[PASS_ALL], na.rm = TRUE),
              #
              N_READS_MUTATED_PTSPEC_ALL_FILTERS = sum(ALT_F[PATIENT_SPECIFIC & PASS_ALL] + ALT_R[PATIENT_SPECIFIC & PASS_ALL]),
              N_READS_MUTATED_NON_PTSPEC_ALL_FILTERS = sum(ALT_F[!PATIENT_SPECIFIC & PASS_ALL] + ALT_R[!PATIENT_SPECIFIC & PASS_ALL]),
              # Number of informative reads (ie total number of reads in a sample, post all filtering)
              N_PTSPEC_INFORMATIVE_READS = sum(REF_F[ALL_IR] + REF_R[ALL_IR] + ALT_F[ALL_IR] + ALT_R[ALL_IR]),
              ALL_INFORMATIVE_READS = sum(REF_F[ALL_IR_CONTROL] + REF_R[ALL_IR_CONTROL] + ALT_F[ALL_IR_CONTROL] + ALT_R[ALL_IR_CONTROL]),
              .groups = "drop") %>%
    left_join(layoutTableTimepoint, by = "SAMPLE_ID") %>%
    left_join(tumourMutationTableSummary, by = "PATIENT") %>%
    full_join(invarScoresTable, by = c('SAMPLE_ID', 'PATIENT')) %>%
    mutate(INITIAL_N_MUTATIONS = ifelse(CASE_OR_CONTROL == 'case', INITIAL_MUTATIONS, 0)) %>%
    mutate(USING_SIZE = ifelse(USING_SIZE, 'WITH_SIZE', 'NO_SIZE'),
           OUTLIER.PASS = ifelse(OUTLIER.PASS, 'PASS', 'FAIL')) %>%
    pivot_wider(names_from = c(USING_SIZE, OUTLIER.PASS),
                names_glue = "DETECTED.{USING_SIZE}.OUTLIER_{OUTLIER.PASS}",
                values_from = DETECTION) %>%
    select(SAMPLE_ID, PATIENT, TIMEPOINT, CASE_OR_CONTROL,
           RAW_N_READS, N_PTSPEC_INFORMATIVE_READS, ALL_INFORMATIVE_READS,
           INITIAL_N_MUTATIONS,
           N_LOCI_MUTATED_PTSPEC, N_LOCI_MUTATED_NON_PTSPEC,
           N_LOCI_MUTATED_PTSPEC_LNP, N_LOCI_MUTATED_NON_PTSPEC_LNP,
           N_LOCI_MUTATED_PTSPEC_ALL_FILTERS, N_LOCI_MUTATED_NON_PTSPEC_ALL_FILTERS,
           N_READS_MUTATED_PTSPEC, N_READS_MUTATED_NON_PTSPEC,
           N_READS_MUTATED_PTSPEC_LNP, N_READS_MUTATED_NON_PTSPEC_LNP,
           N_READS_MUTATED_PTSPEC_ALL_FILTERS, N_READS_MUTATED_NON_PTSPEC_ALL_FILTERS,
           any_of(c('DETECTED.WITH_SIZE.OUTLIER_PASS', 'DETECTED.NO_SIZE.OUTLIER_PASS',
                    'DETECTED.WITH_SIZE.OUTLIER_FAIL', 'DETECTED.NO_SIZE.OUTLIER_FAIL'))) %>%
    arrange(SAMPLE_ID, PATIENT, TIMEPOINT, CASE_OR_CONTROL)

  mutationTracking
}
