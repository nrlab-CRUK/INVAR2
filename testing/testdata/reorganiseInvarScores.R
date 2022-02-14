library(assertthat)
library(tidyverse)

args <- commandArgs(TRUE)

assert_that(length(args) %in% c(2, 4), msg = "Expect two or four arguments: <processed CSV> <layoutFile> [<pool> <barcode>]")

rawFile <- args[1]
layoutFile <- args[2]

toChar <- function(v) { ifelse(v, 'T', 'F') }

layoutTable <-
    suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE)) %>%
    rename_with(str_to_upper) %>%
    rename(POOL = SLX_ID) %>%
    select(POOL, BARCODE, PATIENT, SAMPLE_NAME)

colTypes <- 'diddccllldcldddlciil'

t <- read_csv(rawFile, col_types = colTypes) %>%
    rename_with(str_to_upper)

poolsInData <- 'POOL' %in% colnames(t)
filenameBase <- "REFERENCE_invar_scores"

if (!poolsInData) {
    assert_that(length(args) == 4, msg = "No pool and barcode found in csv. Must be given as arguments.")

    t <- t %>%
        mutate(POOL = args[3], BARCODE = args[4])

    patient <- unique(t$PATIENT)
    assert_that(length(patient) == 1, msg = "More than one patient in the file with explicit pool.")

    filenameBase <- str_c(filenameBase, args[3], args[4], patient, sep = '.')
}

t <- t %>%
    left_join(layoutTable, by = c('POOL', 'BARCODE', 'PATIENT')) %>%
    rename(ALTERNATIVE_LIKELIHOOD = ALT_LIKELIHOOD,
           PATIENT_MUTATION_BELONGS_TO = PT_MUTATION_BELONGS_TO,
           BOTH_STRANDS.PASS = BOTH_STRANDS,
           OUTLIER.PASS = PASS,
           MUTATION_SUM = MUT_SUM) %>%
    select(POOL, BARCODE, SAMPLE_NAME, PATIENT, PATIENT_MUTATION_BELONGS_TO,
           ITERATION, USING_SIZE,
           LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS, CONTAMINATION_RISK.PASS,
           INVAR_SCORE, AF_P, NULL_LIKELIHOOD, ALTERNATIVE_LIKELIHOOD,
           DP, MUTATION_SUM, IMAF, SMOOTH, OUTLIER_SUPPRESSION, MUTANT_READS_PRESENT) %>%
    arrange(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO,
            ITERATION, USING_SIZE,
            LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS)

t %>%
    saveRDS(str_c(filenameBase, ".rds"))

t %>%
    mutate_if(is.double, signif, digits = 6) %>%
    mutate_if(is.logical, toChar) %>%
    write_tsv(str_c(filenameBase, ".tsv"))

t %>%
    filter(ITERATION == 1) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    mutate_if(is.logical, toChar) %>%
    write_tsv(str_c(filenameBase, ".trimmed.tsv"))

