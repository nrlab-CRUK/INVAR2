# This is the code to convert the source file from Emma's run,
# output_size/SLX-19721_SXTLI001.FS_2.sort.filtered.bam.inserts_for_annotation.csv
# into the format the pipeline's code produces.

library(tidyverse)

args <- commandArgs(TRUE)

original <- read_csv(args[1], col_names = c('SNV_BASE', 'SNV', 'STATUS', 'SIZE'), col_types = 'ccci')

fixed <- original %>%
    separate(SNV, into = c('CHROM', 'POS', 'REF', 'ALT'), sep = ':', remove = TRUE) %>%
    mutate(POS = as.integer(POS)) %>%
    select(CHROM, POS, REF, ALT, SNV_BASE, SIZE)

write_tsv(fixed, file = args[2])

