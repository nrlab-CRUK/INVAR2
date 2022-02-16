library(assertthat)
library(tidyverse)

toChar <- function(l)
{
    ifelse(l, 'T', 'F')
}

args <- commandArgs(TRUE)

invisible(assert_that(length(args) == 2, msg = "Expected exactly two arguments: the source table and the output file."))
invisible(assert_that(file.exists(args[1]), msg = str_c(args[1], " does not exist.")))

commonColumnNames <- c("PATIENT", "TIMEPOINT", "CASE_OR_CONTROL",
                       "INITIAL_MUTATIONS", "LOCUS_NOISE.PASS",
                       "MUTANTS_PATIENT_SPECIFIC", "MUTANTS_NON_SPECIFIC",
                       "SPECIFIC_OUTLIER.PASS", "SPECIFIC.PASS", "NON_SPECIFIC.PASS",
                       "WITH_SIZE.OUTLIER_PASS", "NO_SIZE.OUTLIER_PASS",
                       "WITH_SIZE.OUTLIER_FAIL", "NO_SIZE.OUTLIER_FAIL")

originalColumnNames <- c("POOL_BARCODE", commonColumnNames)

currentColumnNames <- c("POOL", "BARCODE", commonColumnNames)

read_csv(args[1], col_names = originalColumnNames, col_types = 'cccciiiiiiillll', skip = 1) %>%
    separate("POOL_BARCODE", into = c("POOL", "BARCODE"), sep = '_') %>%
    select(all_of(currentColumnNames)) %>%
    arrange(POOL, BARCODE, PATIENT, TIMEPOINT, CASE_OR_CONTROL) %>%
    mutate_if(is.logical, toChar) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    write_csv(args[2])
