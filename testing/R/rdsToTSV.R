# Converts an RDS file into a TSV. Has some testing to identify what type of table
# it is likely to be dealing with, and converts accordingly.

library(assertthat)
library(tidyverse)

toChar <- function(x)
{
    ifelse(x, 'T', 'F')
}

mutationTableOrderColumns <-
    c('SAMPLE_ID', 'PATIENT', 'SAMPLE_NAME', 'PATIENT_MUTATION_BELONGS_TO',
      'CHROM', 'POS', 'REF', 'ALT', 'TRINUCLEOTIDE', 'SIZE', 'MUTANT')


args <- commandArgs(TRUE)

invisible(assert_that(length(args) == 2, msg = "Expect exactly two arguments: the RDS file and the output TSV file."))
invisible(assert_that(file.exists(args[1]), msg = str_c(args[1], " does not exist.")))

table <- readRDS(args[1])

if (is.list(table) && !is_tibble(table) && length(table) == 4)
{
    print("Error rate table")

    # Just turn them into one big table

    table <- bind_rows(table)
}
if (all(c('CHROM', 'POS', 'REF', 'ALT', 'DP', 'DP4', 'REF_F', 'ALT_F', 'REF_R', 'ALT_R', 'MQSB') %in% colnames(table)))
{
    print("Mutations table")

    table <- table %>%
        arrange_at(vars(any_of(mutationTableOrderColumns)))
}
if (all(c('N_SAMPLES', 'N_SAMPLES_WITH_SIGNAL') %in% colnames(table)))
{
    print("Locus noise error table")

    table <- table %>%
        arrange(CHROM, POS, TRINUCLEOTIDE)
}

table %>%
    mutate_if(is.logical, toChar) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    write_tsv(args[2])
