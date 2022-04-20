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
    # When reading from an unconverted file, it was:
    # read_csv(tumourMutationsFile, col_types = 'ciccd', show_col_types = FALSE) %>%
    #     select(-contains('uniq'), -any_of('mut')) %>%
    #     rename_with(str_trim) %>%
    #     rename_with(str_to_upper) %>%
    #     rename(CHROM = CHR) %>%
    #     mutate(MUTATION_CLASS = str_c(REF, ALT, sep='/'),
    #            UNIQUE_POS = str_c(CHROM, POS, sep=':'),
    #            UNIQUE_ALT = str_c(UNIQUE_POS, MUTATION_CLASS, sep='_'),
    #            UNIQUE_PATIENT_POS = str_c(PATIENT, UNIQUE_POS, sep='_'))

    read_csv(tumourMutationsFile, col_types = 'ciccd', show_col_types = FALSE) %>%
        mutate(MUTATION_CLASS = str_c(REF, ALT, sep='/'),
               UNIQUE_POS = str_c(CHROM, POS, sep=':'),
               UNIQUE_ALT = str_c(UNIQUE_POS, MUTATION_CLASS, sep='_'),
               UNIQUE_PATIENT_POS = str_c(PATIENT, UNIQUE_POS, sep='_'))
}

loadLayoutTable <- function(layoutFile)
{
    # When reading from an unconverted file, it was:
    # suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE)) %>%
    #     rename_with(str_trim) %>%
    #     rename_with(str_to_upper) %>%
    #     mutate(SAMPLE_ID = str_c(SLX_ID, str_replace(BARCODE, '-', '_'), sep = '_'), .before = "SLX_ID") %>%
    #     select(-SLX_ID, -BARCODE)

    suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, col_types = cols(.default = "c"))) %>%
        mutate(across(INPUT_INTO_LIBRARY_NG, as.double))
}

##
# Typical derived column additions to a mutation table.
##

addMutationTableDerivedColumns <- function(mutationTable)
{
    mutationTable <- mutationTable %>%
        mutate(MUTATION_SUM = ALT_F + ALT_R,
               UNIQUE_POS = str_c(CHROM, POS, sep=':'),
               UNIQUE_ALT = str_c(UNIQUE_POS, str_c(REF, ALT, sep='/'), sep='_'))

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
        select(-any_of(c('MUTATION_SUM', 'PATIENT_SPECIFIC')), -contains('UNIQUE'))
}

##
# Saving functions.
#

# Sorts a mutation table before exporting so the files are in a known order.
# Helps comparison between files.
arrangeMutationTableForExport <- function(mutationTable)
{
    orderByColumns = c('SAMPLE_ID', 'PATIENT', 'PATIENT_MUTATION_BELONGS_TO',
                       'CHROM', 'POS', 'REF', 'ALT', 'TRINUCLEOTIDE', 'SIZE', 'MUTANT')

    mutationTable %>%
        arrange_at(vars(any_of(orderByColumns)))
}

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

# As exportTSV, but write CSV.
exportCSV <- function(t, file)
{
    toChar <- function(x)
    {
        ifelse(x, 'T', 'F')
    }

    t %>%
        mutate_if(is.logical, toChar) %>%
        mutate_if(is.double, signif, digits = 6) %>%
        write_csv(file)

    t
}

# Converts a string so it is safe to use in a file name.
# Replaces all white space with a single underscore, and removes
# all characters that are not word characters (letter, number, underscore).
makeSafeForFileName <- function(string)
{
    string %>%
        str_replace_all("\\s+", "_") %>%
        str_replace_all("[^\\w]+", "")
}
