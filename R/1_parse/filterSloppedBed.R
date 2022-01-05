suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2)
{
    stop("Wrong number of arguments given to filterSloppedBed.R. Expect two.")
}

inFile <- args[1]
outFile <- args[2]

maximumIntervalSize <- 100

read_tsv(inFile,
    col_names = c('CHROM', 'START', 'END', 'REF', 'ALT'),
    col_types = 'ciicc') %>%
filter(END - START <= maximumIntervalSize) %>%
write_tsv(file = outFile, col_names = FALSE)
