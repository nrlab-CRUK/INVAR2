# This is the code to convert the source file from Emma's run,
# output_gz/PARADIGM.f0.9_s2.BQ_20.MQ_40.combined.final.ann.tsv.gz
# into the format expected by the new createMutationsTable.R script.

library(tidyverse)

toChar <- function(l)
{
    ifelse(l, 'T', 'F')
}

TAPAS.tsv_colnames <- c('CHROM', 'POS', 'REF', 'ALT', 'DP', 'DP4','REF_F', 'ALT_F','REF_R', 'ALT_R', 'MQSB', 'POOL', 'BARCODE', 'COSMIC_MUTATIONS', 'COSMIC_SNP', '1KG_AF', 'TRINUCLEOTIDE', 'FILE_NAME')

read_tsv('PARADIGM.f0.9_s2.BQ_20.MQ_40.combined.final.ann.tsv', col_names = TAPAS.tsv_colnames, col_types='ciccicdddddccildcc') %>%
select(-FILE_NAME) %>%
mutate(COSMIC_SNP = toChar(COSMIC_SNP)) %>%
write_tsv('mutation_table.tsv')

