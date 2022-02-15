# Script to convert the four error table off target lists into individual
# mutation table TSV files.
# Doesn't recalculate anything. It brings the column names up to date and
# removes derived columns. It arranges the file by the standard order
# of POOL, BARCODE, CHROM, POS, REF, ALT, TRINULCEOTIDE.


library(assertthat)
library(tidyverse)

loadRData <- function(file)
{
    e <- new.env(parent = emptyenv())
    load(file, envir = e)
    names <- ls(envir = e)
    n <- length(names)
    assert_that(n == 1, msg = str_c(file, " contains ", n, " items: ", paste(names, collapse = ' ')))
    get(names[1], envir = e)
}

toChar <- function(l)
{
    ifelse(l, 'T', 'F')
}

convert <- function(errorRateTable)
{
    as_tibble(errorRateTable) %>%
    rename_all(str_to_upper) %>%
    select(-contains('UNIQ')) %>%
    rename(POOL = SLX, DP_SUM = TOTAL_DP, MUTATION_SUM = MUT_SUM) %>%
    mutate_if(is.logical, toChar) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    arrange(POOL, BARCODE, REF, ALT, TRINUCLEOTIDE)
}

convertAndSave <- function(errorRateTable, tsvFile)
{
    message("Writing ", tsvFile)

    errorRateTable %>%
    convert() %>%
    write_tsv(tsvFile)
}

errorRateMapping <- list(
    PREFILTER = 'pre_filter',
    LOCUS_NOISE = 'locus_noise_filter_only',
    BOTH_STRANDS = 'both_strands_only',
    LOCUS_NOISE.BOTH_STRANDS = 'locus_noise.both_strands')

for (cosmic in c(FALSE, TRUE))
{
    mutationTableFile <- str_c('PARADIGM.f0.9_s2.BQ_20.MQ_40.cosmic_', cosmic, '.error_rates.Rdata')

    message("Reading ", mutationTableFile)

    errorRateList <- loadRData(mutationTableFile)

    cosmicc = ifelse(cosmic, 'cosmic', 'no_cosmic')

    converted <- lapply(errorRateList, convert)
    names(converted) <- names(errorRateMapping)
    saveRDS(converted, str_c('error_rates.off_target.', cosmicc, '.rds'))

    write_tsv(bind_rows(converted), str_c('REFERENCE_error_rates.off_target.', cosmicc, '.tsv'))

    write_tsv(converted$PREFILTER, str_c('REFERENCE_error_rates.off_target.', cosmicc, '.prefilter.tsv'))
    write_tsv(converted$LOCUS_NOISE, str_c('REFERENCE_error_rates.off_target.', cosmicc, '.locusnoise.tsv'))
    write_tsv(converted$BOTH_STRANDS, str_c('REFERENCE_error_rates.off_target.', cosmicc, '.bothstrands.tsv'))
    write_tsv(converted$LOCUS_NOISE.BOTH_STRANDS, str_c('REFERENCE_error_rates.off_target.', cosmicc, '.locusnoise_bothstrands.tsv'))
}

