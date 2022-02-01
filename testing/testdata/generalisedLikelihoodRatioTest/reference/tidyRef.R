library(tidyverse)

toChar <- function(v) { ifelse(v, 'T', 'F') }

read_csv('invarscore.SLX-19721_SXTLI001_PARA_002.csv', col_types = 'dldddccllldldddlii') %>%
    rename_with(str_to_upper) %>%
    mutate(POOL = 'SLX-19721', BARCODE = 'SXTLI001', SAMPLE_NAME = "EXP3079_PARA_002_EOT") %>%
    rename(ALTERNATIVE_LIKELIHOOD = ALT_LIKELIHOOD) %>%
    select(POOL, BARCODE, SAMPLE_NAME, PATIENT, PATIENT_MUTATION_BELONGS_TO,
           LOCUS_NOISE.PASS, BOTH_STRANDS, OUTLIER.PASS, CONTAMINATION_RISK.PASS,
           ITERATION, USING_SIZE,
           INVAR_SCORE, AF_P, NULL_LIKELIHOOD, ALTERNATIVE_LIKELIHOOD,
           DP, MUTATION_SUM, IMAF, SMOOTH, OUTLIER_SUPPRESSION, MUTANT_READS_PRESENT) %>%
    arrange(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO,
            LOCUS_NOISE.PASS, BOTH_STRANDS, OUTLIER.PASS, ITERATION, USING_SIZE) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    mutate_if(is.logical, toChar) %>%
    write_tsv('REFERENCE_invar_scores.SLX-19721_SXTLI001.trimmed.tsv')

