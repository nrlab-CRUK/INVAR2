suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
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
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by part one of the pipeline",
                    default=defaultMarker),
        make_option(c("--pool"), type="character", metavar="string",
                    dest="POOL", help="The pool (SLX) identifier for the BAM/sizes file",
                    default=defaultMarker),
        make_option(c("--barcode"), type="character", metavar="string",
                    dest="BARCODE", help="The barcode for the BAM/sizes file",
                    default=defaultMarker),
        make_option(c("--outlier-suppression"), type="double", metavar="number",
                    dest="OUTLIER_SUPPRESSION", help="The outlier suppression threshold",
                    default=0.05))

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
    list(
        #MUTATIONS_TABLE_FILE = 'on_target/mutation_table.with_sizes.SLX-19721_SXTLI001.rds',
        MUTATIONS_TABLE_FILE = str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/markOutliers/source/combined.polished.size_ann.SLX-19721_SXTLI001.rds'),
        POOL = 'SLX-19721',
        BARCODE = 'SXTLI001',
        OUTLIER_SUPPRESSION = 0.05,
        SAMPLING_SEED = 1024L
    )
}


##
# Calculation functions.
#

##
# From detection_functions.R
#

## Estimate p using the dervied EM algorithm.
# M = MUTANT, R = DP, AF = TUMOUR_AF, e = BACKGROUND_AF
estimate_p_EM <- function(M, R, AF, e, initial_p = 0.01, iterations = 200)
{
    g = AF*(1-e) + (1-AF)*e
    p <- initial_p
    for (i in 1:iterations)
    {
        ## Expectation step
        Z_0 <- (1-g)*p/((1-g)*p + (1-e)*(1-p))
        Z_1 <- g*p/(g*p + e*(1-p))

        ## Maximization step
        p <- sum(M*Z_1 + (R-M)*Z_0) / sum(R)
    }

    p
}


##
# From TAPAS_functions.R
#

# The original version of this method worked on files limited to a single sample.
# This version works on all the samples belonging to a pool + barcode pair,
# performing the sample calculations per grouping of sample and mutation patient.

repolish <- function(mutationTable, outlierSuppressionThreshold)
{
    assert_that(is.number(outlierSuppressionThreshold), msg = "outlierSuppressionThreshold must be a number")

    groupCols <- c('SAMPLE_NAME', 'PATIENT_MUTATION_BELONGS_TO')

    # do not include loci with AF>0.25 or with >10 mutant reads as we are trying to detect MRD.
    # If there are loci with high AF, then there will also be some loci with lower AF (given sufficient loci).
    # If error classes have zero error rate, increase error rate to 1 / BACKGROUND_DP
    # Original code provided a list of 1 for the 'R' parameter. This corresponds to the
    # DP column, so set that on a derived table
    # If there are no mutant reads, AF (p) estimate has to be zero
    # calculate p-value threshold with a bonferroni correction for the number of loci tested

    estimateP <-function(MUTANT, DP, AF, TUMOUR_AF, BACKGROUND_AF)
    {
        max(estimate_p_EM(MUTANT, DP, TUMOUR_AF, BACKGROUND_AF), weighted.mean(AF, TUMOUR_AF))
    }

    mutationTable.forPEstimate <- mutationTable %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS & AF < 0.25 & MUTATION_SUM < 10 & TUMOUR_AF > 0) %>%
        mutate(BACKGROUND_AF = ifelse(BACKGROUND_AF == 0, 1 / BACKGROUND_DP, BACKGROUND_AF),
               DP = 1) %>%
        group_by_at(all_of(groupCols)) %>%
        summarise(P_THRESHOLD = outlierSuppressionThreshold / n_distinct(UNIQUE_POS),
                  P_ESTIMATE = ifelse(!any(MUTANT), 0, estimateP(MUTANT, DP, AF, TUMOUR_AF, BACKGROUND_AF)),
                  .groups = "drop")

    # now take the binom probability of seeing N mutant reads given p (ctDNA estimate in sample)
    # samples pass the filter if they have p-values greater than that defined here.
    # (Passing filter means flagging the mutation as NORMAL, i.e. not an outlier.)

    binom <- function(x, n, p)
    {
        ifelse(x <= 0, 1, binom.test(x, n, p, alternative = "greater")$p.value)
    }

    mutationTable.withTest <- mutationTable %>%
        left_join(mutationTable.forPEstimate, by = groupCols) %>%
        rowwise() %>%
        mutate(BINOMIAL_PROB = binom(x = MUTATION_SUM, n = DP, p = P_ESTIMATE)) %>%
        ungroup() %>% # Removes the "rowwise" grouping
        mutate(OUTLIER.PASS = BINOMIAL_PROB > P_THRESHOLD) %>%
        select(-P_THRESHOLD, -P_ESTIMATE, -BINOMIAL_PROB)

    mutationTable.withTest
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    mutationTable <- mutationTable %>%
        repolish(outlierSuppressionThreshold = scriptArgs$OUTLIER_SUPPRESSION)

    mutationTable %>%
        removeMutationTableDerivedColumns() %>%
        arrangeMutationTableForExport() %>%
        saveRDSandTSV(str_c('mutation_table.outliersuppressed.', scriptArgs$POOL, '_', scriptArgs$BARCODE, '.rds'))
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
