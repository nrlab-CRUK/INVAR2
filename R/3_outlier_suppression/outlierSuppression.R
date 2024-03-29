suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))

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
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) written by sizeAnnotation.R",
                    default=defaultMarker),
        make_option(c("--sample"), type="character", metavar="string",
                    dest="SAMPLE_ID", help="The sample identifier for the mutations table file",
                    default=defaultMarker),
        make_option(c("--patient"), type="character", metavar="string",
                    dest="PATIENT_ID", help="The patient identifier for the mutations table file",
                    default=defaultMarker),
        make_option(c("--outlier-suppression"), type="double", metavar="number",
                    dest="OUTLIER_SUPPRESSION", help="The outlier suppression threshold",
                    default=0.05),
        make_option(c("--allele-frequency-threshold"), type="double", metavar="num",
                    dest="ALLELE_FREQUENCY_THRESHOLD", help="Maximum allele frequency value for acceptable samples",
                    default=0.01),
        make_option(c("--maximum-mutant-reads"), type="integer", metavar="int",
                    dest="MAXIMUM_MUTANT_READS", help="Maximum number of reads acceptable for trying to detect MRD",
                    default=10L))

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


##
# From TAPAS_functions.R
#

repolish <- function(mutationsTable, outlierSuppressionThreshold, alleleFrequencyThreshold, maximumMutantReads)
{
    assert_that(is.number(outlierSuppressionThreshold), msg = "outlierSuppressionThreshold must be a number")

    # do not include loci with AF>alleleFrequencyThreshold or with >maximumMutantReads mutant reads as we are trying to detect MRD.
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

    mutationsTable.forPEstimate <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS & AF <= alleleFrequencyThreshold &
                 MUTATED_READS_PER_LOCI <= maximumMutantReads & TUMOUR_AF > 0) %>%
        mutate(BACKGROUND_AF = ifelse(BACKGROUND_AF == 0, 1 / BACKGROUND_DP, BACKGROUND_AF),
               DP = 1) %>%
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

    mutationsTable.withTest <- mutationsTable %>%
        left_join(mutationsTable.forPEstimate, by = character()) %>%
        rowwise() %>%
        mutate(BINOMIAL_PROB = binom(x = MUTATED_READS_PER_LOCI, n = DP, p = P_ESTIMATE)) %>%
        ungroup() %>% # Removes the "rowwise" grouping
        mutate(OUTLIER.PASS = BINOMIAL_PROB > P_THRESHOLD) %>%
        select(-P_THRESHOLD, -BINOMIAL_PROB) # removing this: 10/05 -P_ESTIMATE

    mutationsTable.withTest
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    # Expect the combination of pool, barcode and patient mutation belongs to to be
    # unique in this file. If the file is empty, that's fine too.

    if (nrow(mutationsTable) == 0)
    {
        warning("There are no rows in ", basename(scriptArgs$MUTATIONS_TABLE_FILE))

        # The table won't have a column for OUTLIER.PASS
        # It should be added before saving.
        mutationsTable <- mutationsTable %>%
            mutate(OUTLIER.PASS = logical())
    }
    else
    {
        mutationsFileCheck <- mutationsTable %>%
            distinct(SAMPLE_ID, PATIENT_MUTATION_BELONGS_TO) %>%
            mutate(MATCHING = SAMPLE_ID == scriptArgs$SAMPLE_ID & PATIENT_MUTATION_BELONGS_TO == scriptArgs$PATIENT_ID)

        assert_that(nrow(mutationsFileCheck) == 1,
                    msg = str_c("Do not have single pool + barcode and patient (mutation belongs to) in ", scriptArgs$MUTATIONS_TABLE_FILE))

        assert_that(all(mutationsFileCheck$MATCHING),
                    msg = str_c("Mutations in", scriptArgs$MUTATIONS_TABLE_FILE, "do not belong to ",
                                scriptArgs$SAMPLE_ID, "and patient (mutation belongs to)",
                                scriptArgs$PATIENT_ID, sep = " "))

        mutationsTable <- mutationsTable %>%
            repolish(outlierSuppressionThreshold = scriptArgs$OUTLIER_SUPPRESSION,
                     alleleFrequencyThreshold = scriptArgs$ALLELE_FREQUENCY_THRESHOLD,
                     maximumMutantReads = scriptArgs$MAXIMUM_MUTANT_READS)

       if (any(is.na(mutationsTable$P_ESTIMATE)))
       {
           stop("Have NA introduced in P_ESTIMATE for ", scriptArgs$SAMPLE_ID, " ", scriptArgs$PATIENT_ID)
       }
    }

    outputName <- str_c('mutation_table.outliersuppressed',
                        makeSafeForFileName(scriptArgs$SAMPLE_ID),
                        makeSafeForFileName(scriptArgs$PATIENT_ID),
                        'rds', sep = '.')

    mutationsTable %>%
        removeMutationTableDerivedColumns() %>%
        arrangeMutationTableForExport() %>%
        saveRDS(outputName)
}

# Launch it.

invisible(main(parseOptions()))
