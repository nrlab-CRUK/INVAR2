suppressWarnings(require(dplyr))
suppressWarnings(require(readr))
suppressWarnings(require(stringr))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5)
{
    stop("Wrong number of arguments given to formatAndAnnotateMutations.R. Expect five.")
}

inFile <- args[1]
outFile <- args[2]
minDP <- as.integer(args[3])
slx <- args[4]
barcode <- args[5]

table <-
    read_tsv(inFile,
        col_names = c('CHROM', 'POS', 'REF', 'ALT', 'DP', 'DP4', 'ADF', 'ADR', 'MQSB'),
        col_types = 'ciccicccc')

converted <- table %>%
    filter(DP >= minDP & MQSB != '.') %>%
    mutate(REF_F = as.integer(ifelse(ALT == '.', ADF, str_split_fixed(ADF, ',', Inf)[,1])),
           ALT_F = as.integer(ifelse(ALT == '.', 0, str_split_fixed(ADF, ',', Inf)[,2])),
           REF_R = as.integer(ifelse(ALT == '.', ADR, str_split_fixed(ADR, ',', Inf)[,1])),
           ALT_R = as.integer(ifelse(ALT == '.', 0, str_split_fixed(ADR, ',', Inf)[,2])),
           .before = 'MQSB') %>%
    select(-ADF, -ADR) %>%
    mutate(MQSB = as.integer(MQSB)) %>%
    mutate(SLX = slx, BARCODE = barcode)

write_tsv(converted, outFile)
