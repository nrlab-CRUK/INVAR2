suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2)
{
    stop("Wrong number of arguments given to patientListCsvToBed.R. Expect two.")
}

csvFile <- args[1]
bedFile <- args[2]

read_csv(csvFile, col_names = TRUE, show_col_types = FALSE) %>%
select(CHROM = chr, POS = pos, REF, ALT) %>%
mutate(POS = as.integer(POS)) %>%
mutate(START = POS - 1, .before = POS) %>%
write_tsv(file = bedFile, col_names = FALSE)
