suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))
source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/detectionFunctions.R'))


##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"

    options_list <- list(
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by outlierSuppression.R",
                    default=defaultMarker),
        make_option(c("--size-characterisation"), type="character", metavar="file",
                    dest="SIZE_CHARACTERISATION_FILE", help="The size characterisation summary file (RDS)",
                    default=defaultMarker),
        make_option(c("--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."),
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
        SIZE_CHARACTERISATION_FILE = str_c(base, '/size_characterisation.rds'),
        BLOODSPOT = FALSE,
        THREADS = 1
    )
}


##
# From TAPAS_functions.R, originally calculate_IMAFv2.
#
calculateIMAFv2 <- function(mutationsTable, bloodspot)
{
    assert_that(is.logical(bloodspot), msg = "bloodspot argument must be a logical")

    mutationsTable.flat <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS, BOTH_STRANDS) %>%
        select(-SIZE, -MUTANT) %>%
        distinct()

    # When it's blood spot data, no need to use outlier suppression as depth is either 1 or 0
    # When not blood spot, also filter on outlier suppression (not an outlier).
    if (!bloodspot)
    {
        mutationsTable.flat <- mutationsTable.flat %>%
            filter(OUTLIER.PASS)
    }

    summary <- mutationsTable.flat %>%
        group_by(MUTATION_CLASS, TRINUCLEOTIDE) %>%
        summarise(TOTAL_DP = sum(REF_F + REF_R + ALT_F + ALT_R),
                  MUTATION_SUM = sum(MUTATION_SUM),
                  MEAN_AF = weighted.mean(AF, DP),
                  BACKGROUND_AF.TRINUCLEOTIDE = first(BACKGROUND_AF),
                  .groups = "drop") %>%
        mutate(MEAN_AF.BS_TRINUCLEOTIDE = pmax(0, MEAN_AF - BACKGROUND_AF.TRINUCLEOTIDE)) %>%
        summarise(IMAFV2 = signif(weighted.mean(MEAN_AF.BS_TRINUCLEOTIDE, TOTAL_DP), 4),
                  .groups = "drop")

    summary$IMAFV2
}

calculateIMAFv2ForAll <- function(mutationsTable, bloodspot)
{
    assert_that(is.logical(bloodspot), msg = "bloodspot argument must be a logical")

    mutationsTable.flat <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS, BOTH_STRANDS) %>%
        select(-SIZE, -MUTANT) %>%
        distinct()

    # When it's blood spot data, no need to use outlier suppression as depth is either 1 or 0
    # When not blood spot, also filter on outlier suppression (not an outlier).
    if (!bloodspot)
    {
        mutationsTable.flat <- mutationsTable.flat %>%
            filter(OUTLIER.PASS)
    }

    summary <- mutationsTable.flat %>%
        group_by(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO, MUTATION_CLASS, TRINUCLEOTIDE) %>%
        summarise(TOTAL_DP = sum(REF_F + REF_R + ALT_F + ALT_R),
                  MUTATION_SUM = sum(MUTATION_SUM),
                  MEAN_AF = weighted.mean(AF, DP),
                  BACKGROUND_AF.TRINUCLEOTIDE = first(BACKGROUND_AF),
                  .groups = "drop") %>%
        group_by(POOL, BARCODE, SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO) %>%
        mutate(MEAN_AF.BS_TRINUCLEOTIDE = pmax(0, MEAN_AF - BACKGROUND_AF.TRINUCLEOTIDE)) %>%
        summarise(IMAFV2 = signif(weighted.mean(MEAN_AF.BS_TRINUCLEOTIDE, TOTAL_DP), 4),
                  .groups = "drop")

    summary
}

## the size table is an aggregate of all the samples, so
## when assessing the detection of a particular sample, exclude itself
## from the aggregated df so that it is not circular
leaveOneOutFilter <- function(mutationsTable, sizeTable)
{
    if (any(mutationsTable$MUTANT))
    {
        readsToSubtract <- mutationsTable %>%
            filter(MUTANT) %>%
            group_by(SIZE) %>%
            summarise(TO_SUBTRACT = n())

        ## do not subtract mutant reads from the wild-type bins
        ## update columns after subtraction

        sizeTable <- sizeTable %>%
            left_join(readsToSubtract, by = "SIZE") %>%
            mutate(TO_SUBTRACT = ifelse(!MUTANT | is.na(TO_SUBTRACT), 0, TO_SUBTRACT)) %>%
            mutate(COUNT = pmax(COUNT - TO_SUBTRACT, 0),
                   TOTAL = pmax(TOTAL - TO_SUBTRACT, 0),
                   PROPORTION = COUNT / TOTAL,
                   TOTAL_MUTANT_READS = sum(ifelse(MUTANT, COUNT, 0))) %>%
            mutate(TOTAL = ifelse(MUTANT, TOTAL_MUTANT_READS, TOTAL)) %>%
            select(-TO_SUBTRACT, -TOTAL_MUTANT_READS)
    }

    sizeTable %>%
        arrange(SIZE)
}


##
# Do the main processing for a single sample and patient mutation pair and
# outlier filter value.
#

doMain <- function(criteria, scriptArgs, mutationsTable, sizeTable)
{
    joinCriteria <- criteria %>%
        select(SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO)

    # Filter down the tables to only be relevant for the criteria.

    mutationsTable <- mutationsTable %>% inner_join(joinCriteria, by = colnames(joinCriteria))

    # Calculate IMAFv2. Now have rows only for one sample and patient mutation
    # is for, so can use a slightly simpler function that returns a single value.

    IMAFv2 <- calculateIMAFv2(mutationsTable, bloodspot = scriptArgs$BLOODSPOT)

    # Apply standard filters to mutation table

    mutationsTable <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS)

    if (criteria$OUTLIER.PASS)
    {
        mutationsTable <- mutationsTable %>%
            filter(OUTLIER.PASS)
    }

    sizeTable <- leaveOneOutFilter(mutationsTable, sizeTable)

    mutationsTable
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    sizeTable <-
        readRDS(scriptArgs$SIZE_CHARACTERISATION_FILE) %>%
        filter(PATIENT_SPECIFIC) %>%
        select(MUTANT, SIZE, COUNT, TOTAL, PROPORTION)

    # Create a table of all samples and variable filter, then turn this into
    # a list of single row tibbles. mclapply can then be used to work on every
    # sample + filter combination in parallel.

    allSamplesAndFilterCombinations <- mutationsTable %>%
        distinct(SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO) %>%
        crossing(OUTLIER.PASS = c(TRUE, FALSE)) %>%
        #filter(PATIENT_MUTATION_BELONGS_TO == 'PARA_028')
        filter(PATIENT_MUTATION_BELONGS_TO == 'PARA_002')

    slicer <- function(n, table) { slice(table, n) }

    allSamplesAndFilterCombinationList <-
        lapply(1:nrow(allSamplesAndFilterCombinations), slicer, allSamplesAndFilterCombinations)

    mclapply(allSamplesAndFilterCombinationList, doMain,
             scriptArgs, mutationsTable, sizeTable,
             mc.cores = scriptArgs$THREADS)

    #IMAFv2
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}