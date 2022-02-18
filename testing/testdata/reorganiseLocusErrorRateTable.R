# Script to convert an error rate data frame into the current format.
# The values are organised by CHROM, POS, TRINUCLEOTIDE.


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

desiredOrder = c('CHROM', 'POS', 'MUTATION_CLASS', 'TRINUCLEOTIDE', 'PATIENT_MUTATION_BELONGS_TO', 'COSMIC',
                 'BACKGROUND_AF', 'MUTATION_SUM', 'DP_SUM', 'N_SAMPLES', 'N_SAMPLES_WITH_SIGNAL', 'LOCUS_NOISE.PASS')

args <- parseArgs(commandArgs(TRUE))

invisible(assert_that(args$EXTENSION %in% names(loadingFunctions), msg = str_c("Unsupported file type: ", args$EXTENSION)))

message("Reading ", args$SOURCE)

table <- loadingFunctions[[args$EXTENSION]](args$SOURCE)

t <- as_tibble(table) %>%
    rename_all(str_to_upper) %>%
    rename(MUTATION_SUM = MUT_SUM, DP_SUM = DP) %>%
    separate(UNIQ_POS, sep = ':', into = c('CHROM', 'POS'), remove = TRUE) %>%
    mutate(POS = as.integer(POS))

if ('MUT_CLASS' %in% colnames(t)) {
    t <- t %>%
        rename(MUTATION_CLASS = MUT_CLASS,
               PATIENT_MUTATION_BELONGS_TO = PT_MUTATION_BELONGS_TO)
}

t <- t %>%
    select(any_of(desiredOrder)) %>%
    arrange(CHROM, POS, TRINUCLEOTIDE)

message("Writing ", args$CONVERTED)

if (args$CONVERTED_EXTENSION == 'rds') {
    saveRDS(t, args$CONVERTED)
} else if (args$CONVERTED_EXTENSION == 'tsv') {
    t %>%
        mutate_if(is.logical, toChar) %>%
        mutate_if(is.double, signif, digits = 6) %>%
        write_tsv(file = args$CONVERTED, progress = TRUE)
}

