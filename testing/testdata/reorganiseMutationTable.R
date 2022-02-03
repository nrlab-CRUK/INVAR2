# Script to convert an original INVAR 1 file containing the mutations table
# into the format created by the INVAR 2 pipeline.
# Doesn't recalculate anything. It brings the column names up to date and
# removes derived columns. It arranges the file by the standard order
# of POOL, BARCODE, CHROM, POS, REF, ALT, TRINULCEOTIDE.


library(assertthat)
library(tidyverse)

loadRData <- function(file)
{
    e <- new.env(parent = emptyenv())
    load(file, envir = e)
    names <- ls(envir = e)
    n <- length(names)
    assert_that(n == 1, msg = str_c(file, " contains ", n, " items: ", paste(names, collapse = ' ')))
    get(names[1], envir = e)
}

toChar <- function(l)
{
    ifelse(l, 'T', 'F')
}

loadingFunctions <- list(rdata = loadRData, rds = readRDS, tsv = read_tsv, csv = read_csv)

parseArgs <- function(args)
{
    assert_that(length(args) >= 1, msg = "Need the name of the mutations table file to convert.")

    filenameParts = (basename(args[1]) %>% str_split('\\.'))[[1]]
    name = str_c(filenameParts[1:length(filenameParts) - 1], collapse='.')
    extension = str_to_lower(tail(filenameParts, 1))

    if (length(args) >= 2)
    {
        convertedFile = args[2]
    }
    else
    {
        dir = dirname(args[1])
        convertedFile = str_c(dir, str_c(name, '.new.tsv'), sep='/')
    }

    filenameParts = (basename(convertedFile) %>% str_split('\\.'))[[1]]
    convertedExtension = str_to_lower(tail(filenameParts, 1))

    list(SOURCE = args[1], NAME = name, EXTENSION = extension, CONVERTED = convertedFile, CONVERTED_EXTENSION = convertedExtension)
}

desiredOrder = c('CHROM', 'POS', 'REF', 'ALT', 'DP', 'DP4', 'REF_F', 'ALT_F', 'REF_R', 'ALT_R', 'MQSB',
                 'POOL', 'BARCODE', 'COSMIC_MUTATIONS', 'COSMIC_SNP', '1KG_AF', 'TRINUCLEOTIDE',
                 'AF', 'COSMIC', 'SNP', 'ON_TARGET',
                 'STUDY', 'SAMPLE_NAME', 'PATIENT', 'SAMPLE_TYPE', 'CASE_OR_CONTROL', 'INPUT_INTO_LIBRARY_NG',
                 'TUMOUR_AF', 'MUTATION_CLASS', 'PATIENT_MUTATION_BELONGS_TO',
                 'BACKGROUND_MUTATION_SUM', 'BACKGROUND_DP', 'BACKGROUND_AF',
                 "LOCUS_NOISE.PASS", "BOTH_STRANDS", "CONTAMINATION_RISK.PASS",
                 "SIZE", "MUTANT", "OUTLIER.PASS")

orderByColumns = c('POOL', 'BARCODE', 'PATIENT', 'SAMPLE_NAME', 'PATIENT_MUTATION_BELONGS_TO',
                   'CHROM', 'POS', 'REF', 'ALT', 'TRINUCLEOTIDE', 'SIZE', 'MUTANT')

args <- parseArgs(commandArgs(TRUE))

invisible(assert_that(args$EXTENSION %in% names(loadingFunctions), msg = str_c("Unsupported file type: ", args$EXTENSION)))

message("Reading ", args$SOURCE)

table <- loadingFunctions[[args$EXTENSION]](args$SOURCE)

message("Converting table")

t <- as_tibble(table) %>%
    rename_all(str_to_upper) %>%
    select(-contains('UNIQ'), -any_of(c('FILE_NAME', 'SLX_BARCODE'))) %>%
    rename(`1KG_AF` = X1KG_AF, POOL = SLX) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(COSMIC_SNP = as.logical(COSMIC_SNP))

if ('BACKGROUND.MUT_SUM' %in% colnames(t)) {
    t <- t %>%
        rename(BACKGROUND_MUTATION_SUM = BACKGROUND.MUT_SUM,
               BACKGROUND_DP = BACKGROUND.DP)
}

if ('PT_MUTATION_BELONGS_TO' %in% colnames(t)) {
    t <- t %>%
        rename(PATIENT_MUTATION_BELONGS_TO = PT_MUTATION_BELONGS_TO)
}

if ('MUT_CLASS' %in% colnames(t)) {
    t <- t %>%
        rename(MUTATION_CLASS = MUT_CLASS)
}

if ('SAMPLE_NAME' %in% colnames(t)) {
    t <- t %>%
        mutate(SAMPLE_NAME = str_remove(SAMPLE_NAME, ' (.+)$'))
}

if ('MUTANT' %in% colnames(t)) {
    t <- t %>%
        mutate(MUTANT = as.logical(MUTANT))
}

if ('PASS' %in% colnames(t)) {
    t <- t %>%
        rename(OUTLIER.PASS = PASS)
}

# See https://stackoverflow.com/questions/26497751/pass-a-vector-of-variable-names-to-arrange-in-dplyr
t <- t %>%
    select(any_of(desiredOrder)) %>%
    arrange_at(vars(any_of(orderByColumns)))

message("Writing ", args$CONVERTED)

if (args$CONVERTED_EXTENSION == 'rds') {
    saveRDS(t, args$CONVERTED)
} else if (args$CONVERTED_EXTENSION == 'tsv') {
    t %>%
        mutate_if(is.logical, toChar) %>%
        mutate_if(is.double, signif, digits = 6) %>%
        write_tsv(t, file = args$CONVERTED)
}


