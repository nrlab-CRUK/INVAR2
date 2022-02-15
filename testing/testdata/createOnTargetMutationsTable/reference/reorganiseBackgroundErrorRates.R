library(assertthat)
library(tidyverse)

toChar <- function(v) { ifelse(v, 'T', 'F') }

args <- commandArgs(TRUE)

t <- readRDS(args[1]) %>%
    rename_with(str_to_upper) %>%
    rename(ERROR_RATE_TYPE = DATA,
           MUTATION_CLASS = MUT_CLASS,
           MUTATION_SUM_TOTAL = BACKGROUND.MUT_SUM,
           TRINUCLEOTIDE_DEPTH = BACKGROUND.DP) %>%
    mutate(ERROR_RATE_TYPE = ifelse(ERROR_RATE_TYPE == 'one_read', 'prefilter', ERROR_RATE_TYPE),
           ERROR_RATE_TYPE = ifelse(ERROR_RATE_TYPE == 'both_reads', 'both_strands', ERROR_RATE_TYPE),
           ERROR_RATE_TYPE = ifelse(ERROR_RATE_TYPE == 'locus_noise.both_reads', 'locus_noise.both_strands', ERROR_RATE_TYPE)) %>%
    arrange(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
    select(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE,
           MUTATION_SUM_TOTAL, TRINUCLEOTIDE_DEPTH, BACKGROUND_AF)

saveRDS(t, 'background_error_rates.rds')

t %>%
    mutate_if(is.logical, toChar) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    write_tsv('REFERENCE_background_error_rates.tsv')

