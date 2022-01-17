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
    extension = str_to_lower(tail(filenameParts, 1))

    if (length(args) >= 2)
    {
        convertedFile = args[2]
    }
    else
    {
        dir = dirname(args[1])
        name = str_c(filenameParts[1:length(filenameParts) - 1], collapse='.')
        convertedFile = str_c(dir, str_c(name, '.new.tsv'), sep='/')
    }

    list(SOURCE = args[1], EXTENSION = extension, CONVERTED = convertedFile)
}

args <- parseArgs(commandArgs(TRUE))

invisible(assert_that(args$EXTENSION %in% names(loadingFunctions), msg = str_c("Unsupported file type: ", args$EXTENSION)))

message("Reading ", args$SOURCE)

table <- loadingFunctions[[args$EXTENSION]](args$SOURCE)

t <- as_tibble(table) %>%
    rename_all(str_to_upper) %>%
    select(-contains('UNIQ'), -any_of(c('FILE_NAME', 'SLX_BARCODE'))) %>%
    rename(`1KG_AF` = X1KG_AF, POOL = SLX) %>%
    mutate(COSMIC_SNP = as.logical(COSMIC_SNP)) %>%
    mutate_if(is.logical, toChar) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    arrange(POOL, BARCODE, CHROM, POS, REF, ALT, TRINUCLEOTIDE)

message("Writing ", args$CONVERTED)

write_tsv(t, file = args$CONVERTED)

