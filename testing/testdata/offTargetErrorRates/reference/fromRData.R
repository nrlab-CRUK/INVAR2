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

convertAndSave <- function(errorRateTable, tsvFile)
{
    message("Writing ", tsvFile)

    as_tibble(errorRateTable) %>%
    rename_all(str_to_upper) %>%
    select(-contains('UNIQ'), -any_of(c('MUT_SUM'))) %>%
    rename(POOL = SLX, DP_SUM = TOTAL_DP) %>%
    mutate_if(is.logical, toChar) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    arrange(POOL, BARCODE, REF, ALT, TRINUCLEOTIDE) %>%
    write_tsv(tsvFile)
}

errorRateMapping <- list(
    ONE_READ = 'pre_filter',
    LOCUS_NOISE = 'locus_noise_filter_only',
    BOTH_READS = 'both_strands_only',
    LOCUS_NOISE.BOTH_READS = 'locus_noise.both_strands')

for (cosmic in c(FALSE, TRUE))
{
    mutationTableFile <- str_c('PARADIGM.f0.9_s2.BQ_20.MQ_40.cosmic_', cosmic, '.error_rates.Rdata')

    message("Reading ", mutationTableFile)

    errorRateList <- loadRData(mutationTableFile)

    cosmicc = ifelse(cosmic, 'cosmic', 'no_cosmic')

    saveRDS(errorRateList, str_c('mutation_table.off_target.', cosmicc, '.rds'))

    convertAndSave(errorRateList$pre_filter, str_c('REFERENCE_mutation_table.off_target.', cosmicc, '.oneread.tsv'))
    convertAndSave(errorRateList$locus_noise_filter_only, str_c('REFERENCE_mutation_table.off_target.', cosmicc, '.locusnoise.tsv'))
    convertAndSave(errorRateList$both_strands_only, str_c('REFERENCE_mutation_table.off_target.', cosmicc, '.bothreads.tsv'))
    convertAndSave(errorRateList$locus_noise.both_strands, str_c('REFERENCE_mutation_table.off_target.', cosmicc, '.locusnoise_bothreads.tsv'))
}

