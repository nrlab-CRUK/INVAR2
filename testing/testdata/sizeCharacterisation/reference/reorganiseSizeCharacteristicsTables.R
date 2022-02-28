# Script to convert the size characteristics tables from existing data
# into the form output by sizeClassification.R

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

parseArg <- function(arg)
{
    filenameParts = (basename(arg) %>% str_split('\\.'))[[1]]
    name = str_c(filenameParts[1:length(filenameParts) - 1], collapse='.')
    extension = str_to_lower(tail(filenameParts, 1))

    convertedFile = str_c(name, '.tsv')

    list(SOURCE = arg, NAME = name, EXTENSION = extension, CONVERTED = convertedFile)
}

desiredOrder = c('SAMPLE_ID', 'PATIENT_MUTATION_BELONGS_TO',
                 'CASE_OR_CONTROL', 'PATIENT_SPECIFIC', 'MUTANT', 'SIZE', 'TOTAL' )

orderByColumns = c('SAMPLE_ID', 'PATIENT_MUTATION_BELONGS_TO', 'MUTANT', 'SIZE')

args <- commandArgs(TRUE)

invisible(assert_that(length(args) == 2, msg = "Expect exactly two files: the all sizes file and the summary file"))

allSizesFile <- parseArg(args[1])

invisible(assert_that(allSizesFile$EXTENSION %in% names(loadingFunctions), msg = str_c("Unsupported file type: ", allSizesFile$EXTENSION)))

message("Reading ", allSizesFile$SOURCE)

table <- loadingFunctions[[allSizesFile$EXTENSION]](allSizesFile$SOURCE)

message("Converting table")

allSizes <- as_tibble(table) %>%
    mutate_if(is.factor, as.character) %>%
    separate(SLX_barcode, into = c('POOL', 'BARCODE'), sep = "_") %>%
    mutate(SAMPLE_ID = str_c(POOL, BARCODE, sep = ":")) %>%
    separate(sample_name, into = c('SAMPLE_NAME', 'PATIENT_MUTATION_BELONGS_TO'), sep = " ") %>%
    mutate(PATIENT_MUTATION_BELONGS_TO = str_remove_all(PATIENT_MUTATION_BELONGS_TO, "[()]")) %>%
    mutate(PATIENT_SPECIFIC = data == 'ptspec') %>%
    rename_all(str_to_upper) %>%
    select(any_of(desiredOrder)) %>%
    arrange_at(vars(any_of(orderByColumns)))

message("Writing size_characterisation.all.rds")

allSizes %>%
    saveRDS('size_characterisation.all.rds')

message("Writing size_characterisation.all.tsv")

allSizes %>%
    mutate_if(is.logical, toChar) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    write_tsv(file = 'REFERENCE_size_characterisation.all.tsv')



summaryFile <- parseArg(args[2])

invisible(assert_that(summaryFile$EXTENSION %in% names(loadingFunctions), msg = str_c("Unsupported file type: ", summaryFile$EXTENSION)))

message("Reading ", summaryFile$SOURCE)

table <- loadingFunctions[[summaryFile$EXTENSION]](summaryFile$SOURCE)

message("Converting summary table")

summary <- as_tibble(table) %>%
    rename_all(str_to_upper) %>%
    mutate_if(is.factor, as.character) %>%
    rename(MUTANT = MUT) %>%
    mutate(PATIENT_SPECIFIC = DATA == 'ptspec') %>%
    select(PATIENT_SPECIFIC, CASE_OR_CONTROL, MUTANT, SIZE, COUNT, TOTAL, PROPORTION) %>%
    arrange(PATIENT_SPECIFIC, CASE_OR_CONTROL, MUTANT, SIZE)

message("Writing size_characterisation.rds")

summary %>%
    saveRDS('size_characterisation.rds')

message("Writing size_characterisation.tsv")

summary %>%
    mutate_if(is.logical, toChar) %>%
    mutate_if(is.double, signif, digits = 6) %>%
    write_tsv(file = 'REFERENCE_size_characterisation.tsv')

