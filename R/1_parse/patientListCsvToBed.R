suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2)
{
    stop("Wrong number of arguments given to patientListCsvToBed.R. Expect two.")
}

csvFile <- args[1]
bedFile <- args[2]

loadTumourMutationsTable(csvFile) %>%
select(CHROM, POS, REF, ALT) %>%
mutate(START = POS - 1, .before = POS) %>%
write_tsv(file = bedFile, col_names = FALSE)
