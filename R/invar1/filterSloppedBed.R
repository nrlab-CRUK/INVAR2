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

table <-
    read_tsv(inFile,
        col_names = c('CHROM', 'START', 'END', 'REF', 'ALT'),
        col_types = 'ciicc')

filtered <- table %>%
    filter(END - START <= maximumIntervalSize)

write_tsv(filtered, outFile, col_names = FALSE)
