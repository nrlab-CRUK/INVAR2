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
                    dest="OUTLIER_SUPPRESSION", help="The outlier suppression setting",
                    default=0.05),
        make_option(c("--sampling-seed"), type="integer", metavar="number",
                    dest="SAMPLING_SEED", help="The seed for random sampling. Only use in testing.",
                    default=NA))

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
        MUTATIONS_TABLE_FILE = str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/outlierSuppression/source/EXP3079_PARA_002_EOT.SLX-19721_SXTLI001.1.ptspec.0.05.combined.polished.size.rds'),
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
estimate_p_EM <- function(mutationTable, initial_p = 0.01, iterations = 200)
{
    M = mutationTable$MUTANT
    R = mutationTable$DP
    AF = mutationTable$TUMOUR_AF
    e = mutationTable$BACKGROUND_AF

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


repolish <- function(mutationTable, outlierSuppressionThreshold)
{
    assert_that(is.number(outlierSuppressionThreshold), msg = "outlierSuppressionThreshold must be a number")

    # do not include loci with AF>0.25 or with >10 mutant reads as we are trying to detect MRD.
    # If there are loci with high AF, then there will also be some loci with lower AF (given sufficient loci).
    # If error classes have zero error rate, increase error rate to 1 / BACKGROUND_DP

    mutationTable.forPEstimate <- mutationTable %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS & AF < 0.25 & MUTATION_SUM < 10 & TUMOUR_AF > 0) %>%
        mutate(BACKGROUND_AF = ifelse(BACKGROUND_AF == 0, 1 / BACKGROUND_DP, BACKGROUND_AF))

    mutationTable.forPEstimate.nUniquePositions <- mutationTable.forPEstimate %>%
        summarise(N_UNIQUE_POS = n_distinct(UNIQUE_POS)) %>%
        as.integer()

    # if no mutant reads, AF estimate has to be zero
    if (!any(mutationTable.forPEstimate$MUTANT))
    {
        p_estimate <- 0
    }
    else
    {
        # Original provided a list of 1 for the 'R' parameter. This corresponds to the
        # DP column, so set that on a derived table

        mutationTable.P <- mutationTable.forPEstimate %>%
            mutate(DP = 1)

        # estimate p for outlier suppression
        p_estimate <- estimate_p_EM(mutationTable.P)
        p_estimate <- max(p_estimate, weighted.mean(mutationTable.P$AF, mutationTable.P$TUMOUR_AF))
    }

    # calculate p-value threshold with a bonferroni correction for the number of loci tested
    p_threshold <- outlierSuppressionThreshold / mutationTable.forPEstimate.nUniquePositions


    binom <- function(x, n, p)
    {
        ifelse(x <= 0, 1, binom.test(x, n, p, alternative = "greater")$p.value)
    }

    # now take the binom probability of seeing N mutant reads given p (ctDNA estimate in sample)
    # samples pass the filter if they have p-values greater than that defined here.
    mutationTable <- mutationTable %>%
        rowwise() %>%
        mutate(BINOMIAL_PROB = binom(x = MUTATION_SUM, n = DP, p = p_estimate)) %>%
        ungroup() %>% # Removes the "rowwise" grouping
        mutate(REPOLISH.PASS = BINOMIAL_PROB > p_threshold) %>%
        select(-BINOMIAL_PROB)

    mutationTable
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns() %>%
        mutate(MUTATION_SUM = ALT_F + ALT_R)

    mutationsTable <- mutationsTable %>%
        repolish(outlierSuppressionThreshold = scriptArgs$OUTLIER_SUPPRESSION)

    mutationsTable %>%
        saveRDSandTSV('mutation_table.repolished.rds')

    ## TESTING

    #aSample <- mutationsTable %>%
    #    filter(PATIENT_SPECIFIC) %>%
    #    distinct(SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO) %>%
    #    slice_head(n = 1)

    # EXP3079_PARA_002_EOT_(PARA_002)

    #mutationsTable2 <- mutationsTable %>%
    #    filter(PATIENT_SPECIFIC &
    #           PATIENT_MUTATION_BELONGS_TO == "PARA_002" &
    #           SAMPLE_NAME == "EXP3079_PARA_002_EOT")

    #repolish(mutationsTable2, scriptArgs$OUTLIER_SUPPRESSION)
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
