# Script to convert an error rate data frame into the current format.
# The values are organised by UNIQUE_POS, TRINUCLEOTIDE.


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
    rename(MUTATION_SUM = MUT_SUM, DP_SUM = DP) %>%
    separate(UNIQ_POS, sep = ':', into = c('CHROM', 'POS'), remove = TRUE) %>%
    mutate(POS = as.integer(POS)) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    arrange(CHROM, POS, TRINUCLEOTIDE) %>%
    select(any_of(c('CHROM', 'POS', 'TRINUCLEOTIDE', 'BACKGROUND_AF', 'MUTATION_SUM', 'DP_SUM', 'N_SAMPLES', 'N_SAMPLES_WITH_SIGNAL')))

message("Writing ", args$CONVERTED)

write_tsv(t, file = args$CONVERTED)

