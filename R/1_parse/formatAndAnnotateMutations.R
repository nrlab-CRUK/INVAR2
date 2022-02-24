suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4)
{
    stop("Wrong number of arguments given to formatAndAnnotateMutations.R. Expect four.")
}

inFile   <- args[1]
outFile  <- args[2]
minDepth <- as.integer(args[3])
sampleId <- args[4]

# Original R code implies in a comment that the filter should also include
# "MQSB != '.'" but it is not present in the actual active code.

read_tsv(inFile,
    col_names = c('CHROM', 'POS', 'REF', 'ALT', 'DP', 'DP4', 'ADF', 'ADR', 'MQSB'),
    col_types = 'ciccicccc') %>%
filter(DP >= minDepth) %>%
mutate(REF_F = as.integer(ifelse(ALT == '.', ADF, str_split_fixed(ADF, ',', Inf)[,1])),
       ALT_F = as.integer(ifelse(ALT == '.', 0, str_split_fixed(ADF, ',', Inf)[,2])),
       REF_R = as.integer(ifelse(ALT == '.', ADR, str_split_fixed(ADR, ',', Inf)[,1])),
       ALT_R = as.integer(ifelse(ALT == '.', 0, str_split_fixed(ADR, ',', Inf)[,2])),
       .before = 'MQSB') %>%
select(-ADF, -ADR) %>%
mutate(MQSB = as.double(MQSB),
       SAMPLE_ID = sampleId) %>%
write_tsv(file = outFile)
