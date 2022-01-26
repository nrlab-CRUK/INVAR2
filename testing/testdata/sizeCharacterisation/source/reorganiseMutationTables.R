# Script to convert an original INVAR 1 file containing the mutations table
# into the format created by the INVAR 2 pipeline.
# Doesn't recalculate anything. It brings the column names up to date and
# removes derived columns. It arranges the file by the standard order
# of POOL, BARCODE, CHROM, POS, REF, ALT, TRINULCEOTIDE.


library(assertthat)
library(parallel)
library(tidyverse)

loadTumourMutationsTable <- function(tumourMutationsFile)
{
    read_csv(tumourMutationsFile, col_types = 'ciccd', show_col_types = FALSE, progress = TRUE) %>%
        rename_with(str_to_upper) %>%
        rename(POOL = SLX_ID) %>%
        select(POOL, BARCODE, CASE_OR_CONTROL)
}

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

parseOne <- function(arg)
{
    filenameParts = (arg %>% str_split('\\.'))[[1]]
    name = str_c(filenameParts[1:length(filenameParts) - 1], collapse='.')
    extension = str_to_lower(tail(filenameParts, 1))

    dir = dirname(arg)
    convertedFile = str_c(dir, str_c(name, '.rds'), sep='/')

    list(SOURCE = arg, NAME = name, EXTENSION = extension, CONVERTED = convertedFile)
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

processOne <- function(arg, tumourMutationsTable)
{
    args <- parseOne(arg)

    message("Reading ", args$SOURCE)

    table <- loadingFunctions[[args$EXTENSION]](args$SOURCE)

    message("Converting table")

    t <- as_tibble(table) %>%
        rename_all(str_to_upper) %>%
        select(-contains('UNIQ'), -any_of(c('FILE_NAME', 'SLX_BARCODE'))) %>%
        rename(`1KG_AF` = X1KG_AF, POOL = SLX) %>%
        mutate_if(is.factor, as.character) %>%
        mutate(COSMIC_SNP = as.logical(COSMIC_SNP)) %>%
        left_join(tumourMutationsTable, by = c('POOL', 'BARCODE'))

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
            mutate(SAMPLE_NAME = str_remove(SAMPLE_NAME, str_c(' (.+)$')))
    }

    if ('MUTANT' %in% colnames(t)) {
        t <- t %>%
            mutate(MUTANT = as.logical(MUTANT))
    }

    if ('PASS' %in% colnames(t)) {
        t <- t %>%
            rename(OUTLIER.PASS = PASS)
    }

    t %>%
        select(any_of(desiredOrder)) %>%
        arrange_at(vars(any_of(orderByColumns)))
}

args <- commandArgs(TRUE)
#args <- c('/home/data/INVAR/source_files/combined.SLX_table_with_controls_031220.csv')
invisible(assert_that(length(args) == 1, msg = "Expect exactly one argument: the path to the layout file (SLX table)."))

tumourMutationsFile <- args[1]
invisible(assert_that(file.exists(tumourMutationsFile), msg = str_c(tumourMutationsFile, " doesn't exist.")))

rdataFiles = list.files(pattern = "\\.Rdata$")
invisible(assert_that(length(rdataFiles) > 0, msg = str_c("No Rdata files found in ", getwd())))

dataFileGrouping <-
    tibble(FILENAME = rdataFiles) %>%
    mutate(POOL_BARCODE = str_match(FILENAME, ".*(SLX-\\d+_\\w+)\\..*")[,2])

tumourMutationsTable <- loadTumourMutationsTable(tumourMutationsFile)

coresToUse <- ceiling(detectCores(logical = FALSE) / 2)
coresToUse <- ifelse(is.na(coresToUse), 1, coresToUse)

message("Processing files using ", coresToUse, " CPUs.")

poolBarcodePairs <- unique(dataFileGrouping$POOL_BARCODE)

for (poolBarcode in poolBarcodePairs) {
    message("Handling files belonging to ", poolBarcode, ". ", match(poolBarcode, poolBarcodePairs), " of ", length(poolBarcodePairs))

    name <- str_c(poolBarcode, '.os.rds')

    if (file.exists(name))
    {
        message(name, " already exists, so skipping.")
    }
    else
    {
        theseFiles <- dataFileGrouping %>%
            filter(POOL_BARCODE == poolBarcode)

        theseTables <- mclapply(theseFiles$FILENAME, processOne, tumourMutationsTable, mc.cores = coresToUse)

        combined <- bind_rows(theseTables)

        message("Writing ", name)
        saveRDS(combined, file = name)
    }
}
