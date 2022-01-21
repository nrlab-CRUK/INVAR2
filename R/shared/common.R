##
# This script is not a stand alone piece of work but contains functions that
# are common across all the invar34 scripts. Some functions would be common to all
# scripts, and others (the loading ones) offer the most basic common part of
# the manipulation of the files. Scripts that require more should extend their
# functionality.
#


##
# Loading functions.
#

loadTumourMutationsTable <- function(tumourMutationsFile)
{
    read_csv(tumourMutationsFile, col_types = 'ciccd', show_col_types = FALSE) %>%
        select(-contains('uniq'), -any_of('mut')) %>%
        rename_with(str_to_upper) %>%
        rename(CHROM = CHR) %>%
        mutate(MUTATION_CLASS = str_c(REF, ALT, sep='/'),
               UNIQUE_POS = str_c(CHROM, POS, sep=':'),
               UNIQUE_ALT = str_c(UNIQUE_POS, MUTATION_CLASS, sep='_'),
               UNIQUE_PATIENT_POS = str_c(PATIENT, UNIQUE_POS, sep='_'))
}

loadLayoutTable <- function(layoutFile)
{
    suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE)) %>%
        rename_with(str_to_upper) %>%
        mutate(POOL_BARCODE = str_c(SLX_ID, str_replace(BARCODE, '-', '_'), sep = '_'))
}

##
# Typical derived column additions to a mutation table.
##

addMutationTableDerivedColumns <- function(mutationTable)
{
    mutationTable <- mutationTable %>%
        mutate(MUTATION_SUM = ALT_F + ALT_R,
               UNIQUE_POS = str_c(CHROM, POS, sep=':'),
               UNIQUE_ALT = str_c(UNIQUE_POS, str_c(REF, ALT, sep='/'), sep='_'),
               POOL_BARCODE = str_c(POOL, BARCODE, sep='_'))

    if (all(c('PATIENT', 'PATIENT_MUTATION_BELONGS_TO') %in% colnames(mutationTable)))
    {
        mutationTable <- mutationTable %>%
            mutate(PATIENT_SPECIFIC = PATIENT == PATIENT_MUTATION_BELONGS_TO)
    }

    mutationTable
}



# Remove columns from the mutation table that can be derived from
# other columns, typically before saving.
removeMutationTableDerivedColumns <- function(mutationTable)
{
    mutationTable %>%
        select(-any_of(c('MUTATION_SUM', 'POOL_BARCODE', 'PATIENT_SPECIFIC')), -contains('UNIQUE'))
}

##
# Saving functions.
#


# Writes the TSV file but, before saving, converts any logical columns to
# simply the characters 'T' or 'F'. Saves having full "TRUE" and "FALSE" values,
# which are excessive as reading the table back correctly interprets 'T' and 'F'.
# Also limits floating point numbers to six significant figures.
exportTSV <- function(t, file)
{
    toChar <- function(x)
    {
        ifelse(x, 'T', 'F')
    }

    t %>%
        mutate_if(is.logical, toChar) %>%
        mutate_if(is.double, signif, digits = 6) %>%
        write_tsv(file)

    t
}

# Save the given table as an RDS file and a TSV
saveRDSandTSV <- function(t, file)
{
    saveRDS(t, file)

    tsv <- str_replace(file, "\\.rds$", ".tsv")

    if (tsv != file)
    {
        exportTSV(t, tsv)
    }

    t
}
