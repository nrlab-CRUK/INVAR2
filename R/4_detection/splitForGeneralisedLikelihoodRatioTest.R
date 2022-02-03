# Split outlier suppressed mutations files by sample.
# Allows the generalised likelihood ratio test to work on a smaller data
# set per job, which on a cluster at least allows more to be done in
# parallel.
# The mutations file coming in is for a pool + barcode pair, not the
# whole thing.

suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))


##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"

    options_list <- list(
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) with outlier suppression flags",
                    default=defaultMarker),
        make_option(c("--pool"), type="character", metavar="string",
                    dest="POOL", help="The pool (SLX) identifier for the mutations table file (one pool + barcode)",
                    default=defaultMarker),
        make_option(c("--barcode"), type="character", metavar="string",
                    dest="BARCODE", help="The barcode for the mutations table file (one pool + barcode)",
                    default=defaultMarker),
        make_option(c("--threads"), type="integer", metavar="integer",
                    dest="THREADS", help="The number of cores to use to process the input files.",
                    default=1))

    opts <- OptionParser(option_list=options_list, usage="%prog [options]") %>%
        parse_args(positional_arguments = TRUE)

    missingRequired <- which(opts$options == defaultMarker)
    if (any(missingRequired))
    {
        missing <- names(opts$options)[missingRequired]
        stop(str_c("ERROR: One or more required arguments have not been provided: ", str_c(missing, collapse = ' ')))
    }

    scriptOptions <- list()
    for (v in names(opts$options))
    {
        scriptOptions[v] = opts$options[v]
    }

    scriptOptions
}

# Test options for my (Rich) local set up in RStudio.

richTestOptions <- function()
{
    base <- str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/generalisedLikelihoodRatioTest/source')

    list(
        MUTATIONS_TABLE_FILE = str_c(base, '/SLX-19721_SXTLI001.os.rds'),
        POOL = 'SLX-19721',
        BARCODE = 'SXTLI001',
        THREADS = 4L
    )
}

##
# A per sample main method. Filter on sample name and save.
#

doMain <- function(patientMutationBelongsTo, scriptArgs, mutationsTable)
{
    mutationsTable <- mutationsTable %>%
        filter(PATIENT_MUTATION_BELONGS_TO == patientMutationBelongsTo)

    fnPatientName <- str_replace_all(patientMutationBelongsTo, "\\s+", "_")
    fnPatientName <- str_replace_all(patientMutationBelongsTo, "[^\\w]+", "")

    filename <- str_c(scriptArgs$POOL, "_", scriptArgs$BARCODE, "_", fnPatientName, ".perpatient.rds")

    saveRDS(mutationsTable, filename)

    tibble(POOL = scriptArgs$POOL, BARCODE = scriptArgs$BARCODE,
           PATIENT_MUTATION_BELONGS_TO = patientMutationBelongsTo, FILE_NAME = filename)
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE)

    patientMutations = mutationsTable %>%
        distinct(PATIENT_MUTATION_BELONGS_TO)

    fileInfoList <-
        mclapply(patientMutations$PATIENT_MUTATION_BELONGS_TO, doMain,
                 scriptArgs, mutationsTable, mc.cores = scriptArgs$THREADS)

    fileInfoTable <- bind_rows(fileInfoList)

    if (any(duplicated(fileInfoTable$FILE_NAME)))
    {
        print(fileInfoTable)
        stop("Some file names have clashed after munging. This is a problem that cannot be fixed without changing the patient names.")
    }

    write_csv(fileInfoTable, "index.csv")
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
