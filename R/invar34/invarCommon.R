##
# Loading functions.
#


# Load the patient specific tumour mutations file

load.tumour.mutations.table <- function(tumourMutationsFile)
{
    read_csv(tumourMutationsFile, col_types = 'ciccd', show_col_types = FALSE) %>%
    select(-contains('uniq'), -mut) %>%
    rename_with(str_to_upper) %>%
    rename(CHROM = CHR) %>%
    mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':'))
}

