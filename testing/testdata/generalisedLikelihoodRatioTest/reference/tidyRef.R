library(tidyverse)

toChar <- function(v) { ifelse(v, 'T', 'F') }

standardise <- function(table) {
    table %>%
        rename_with(str_to_upper) %>%
        rename(INVAR_SCORE = INVAR_SCORE.USING_SIZE) %>%
        rename(ALTERNATIVE_LIKELIHOOD = ALT_LIKELIHOOD,
               PATIENT_MUTATION_BELONGS_TO = PT_MUTATION_BELONGS_TO,
               OUTLIER.PASS = PASS,
               MUTATION_SUM = MUT_SUM) %>%
        mutate(USING_SIZE = row_number() %% 2 == 1) %>%
        select(POOL, BARCODE, SAMPLE_NAME, PATIENT, PATIENT_MUTATION_BELONGS_TO,
               ITERATION, USING_SIZE,
               LOCUS_NOISE.PASS, BOTH_STRANDS, OUTLIER.PASS, CONTAMINATION_RISK.PASS,
               INVAR_SCORE, AF_P, NULL_LIKELIHOOD, ALTERNATIVE_LIKELIHOOD,
               DP, MUTATION_SUM, IMAF, SMOOTH, OUTLIER_SUPPRESSION, MUTANT_READS_PRESENT) %>%
        arrange(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO,
                ITERATION, USING_SIZE,
                LOCUS_NOISE.PASS, BOTH_STRANDS, OUTLIER.PASS) %>%
        mutate_if(is.double, signif, digits = 6) %>%
        mutate_if(is.logical, toChar)
}

colTypes <- 'diddccllldcldddlcii'

read_csv('invarscore.SLX-19721_SXTLI001_PARA_002.csv', col_types = colTypes) %>%
    mutate(POOL = 'SLX-19721', BARCODE = 'SXTLI001', SAMPLE_NAME = "EXP3079_PARA_002_EOT") %>%
    standardise() %>%
    write_tsv('REFERENCE_invar_scores.SLX-19721.SXTLI001.PARA_002.tsv')

p028 <-
    read_csv('invarscore.SLX-19721_SXTLI001_PARA_028.csv', col_types = colTypes) %>%
        mutate(POOL = 'SLX-19721', BARCODE = 'SXTLI001', SAMPLE_NAME = "EXP3079_PARA_002_EOT") %>%
        standardise()

p028 %>%
    write_tsv('REFERENCE_invar_scores.SLX-19721.SXTLI001.PARA_028.tsv')

p028 %>%
    filter(ITERATION == 1) %>%
    write_tsv('REFERENCE_invar_scores.SLX-19721.SXTLI001.PARA_028.trimmed.tsv')

