suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
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
                    dest="MUTATIONS_TABLE_FILE", help="The per-patient mutations table file (RDS) with outlier suppression values",
                    default=defaultMarker),
        make_option(c("--size-characterisation"), type="character", metavar="file",
                    dest="SIZE_CHARACTERISATION_FILE", help="The size characterisation summary file (RDS)",
                    default=defaultMarker),
        make_option(c("--sample"), type="character", metavar="string",
                    dest="SAMPLE_ID", help="The sample identifier for the mutations table file",
                    default=defaultMarker),
        make_option(c("--patient"), type="character", metavar="string",
                    dest="PATIENT_ID", help="The patient identifier for the mutations table file",
                    default=defaultMarker),
        make_option(c("--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."),
        make_option(c("--outlier-suppression"), type="double", metavar="number",
                    dest="OUTLIER_SUPPRESSION", help="The outlier suppression threshold",
                    default=0.05),
        make_option(c("--threads"), type="integer", metavar="integer",
                    dest="THREADS", help="The number of cores to use to process the input files.",
                    default=1L),
        make_option(c("--minimum-fragment-length"), type="integer", metavar="integer",
                    dest="MINIMUM_FRAGMENT_LENGTH", help="The minimum fragment length.",
                    default=60L),
        make_option(c("--maximum-fragment-length"), type="integer", metavar="integer",
                    dest="MAXIMUM_FRAGMENT_LENGTH", help="The maximum fragment length.",
                    default=300L),
        make_option(c("--smoothing"), type="double", metavar="number",
                    dest="SMOOTHING", help="Width of smoothing.",
                    default=0.03),
        make_option(c("--only-weigh-mutants"), action="store_true", default=FALSE,
                    dest="ONLY_WEIGH_MUTANTS", help="Only weigh ctDNA signal based on mutant fragments."),
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


##
# From TAPAS_functions.R, originally calculate_IMAFv2.
#
calculateIMAFv2 <- function(mutationsTable, bloodspot)
{
    assert_that(is.logical(bloodspot), msg = "bloodspot argument must be a logical")

    mutationsTable.flat <- mutationsTable %>%
        filter(LOCUS_NOISE.PASS, BOTH_STRANDS.PASS) %>%
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
        summarise(IMAFV2 = weighted.mean(MEAN_AF.BS_TRINUCLEOTIDE, TOTAL_DP),
                  .groups = "drop")

    summary$IMAFV2
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

likelihoodRatioListToTibble <- function(likelihoodRatio)
{
    as_tibble(likelihoodRatio) %>%
        rename_at(vars(ends_with(".no_size")), ~ str_replace_all(., '\\.no_size', '')) %>%
        rename_with(str_to_upper)
}

calculateLikelihoodRatioForSampleWithSize <- function(mutationsTable, sizeTable,
                                                      minFragmentLength, maxFragmentLength,
                                                      smooth, onlyWeighMutants)
{
    assert_that(is.numeric(minFragmentLength), msg = "minFragmentLength must be an number.")
    assert_that(is.numeric(maxFragmentLength), msg = "maxFragmentLength must be an number.")
    assert_that(minFragmentLength <= maxFragmentLength, msg = "minFragmentLength must be <= maxFragmentLength.")
    assert_that(is.numeric(smooth), msg = "smooth must be an number.")
    assert_that(is.logical(onlyWeighMutants), msg = "onlyWeighMutants must be a logical.")

    # Use depth = 1 for these calculations.

    mutationsTable <- mutationsTable %>%
        mutate(DP = 1)

    # Filter dataframe to include only minFragmentLength:maxFragmentLength fragments - affects the weight normalisation
    sizeTable <- filter(sizeTable, SIZE>=minFragmentLength, SIZE<=maxFragmentLength)

    ## read length probabilites for mutant reads

    mutantSizeTable <- sizeTable %>%
        filter(MUTANT)

    mutantProbabilities <-
        estimate_real_length_probability(fragment_length = mutantSizeTable$SIZE,
                                         counts = mutantSizeTable$COUNT,
                                         min_length = minFragmentLength,
                                         max_length = maxFragmentLength,
                                         bw_adjust = smooth) %>%
        as_tibble() %>%
        rename_with(str_to_upper)

    ## read length probabilties of normal reads

    normalSizeTable <- sizeTable %>%
        filter(!MUTANT)

    normalProbabilities <-
        estimate_real_length_probability(fragment_length = normalSizeTable$SIZE,
                                         counts = normalSizeTable$COUNT,
                                         min_length = minFragmentLength,
                                         max_length = maxFragmentLength,
                                         bw_adjust = smooth) %>%
        as_tibble() %>%
        rename_with(str_to_upper)

    # Add these probabilities to the mutations table as new columns.

    mutationsTable <- mutationsTable %>%
        left_join(normalProbabilities, by = c('SIZE' = 'FRAGMENT_LENGTH')) %>%
        rename(REAL_LENGTH_PROB_NORMAL = PROBABILITY) %>%
        left_join(mutantProbabilities, by = c('SIZE' = 'FRAGMENT_LENGTH')) %>%
        rename(REAL_LENGTH_PROB_MUTANT = PROBABILITY)

    assert_that(!any(is.na(mutationsTable$REAL_LENGTH_PROB_NORMAL)), msg = "NAs in normal real length probability")
    assert_that(!any(is.na(mutationsTable$REAL_LENGTH_PROB_MUTANT)), msg = "NAs in mutant real length probability")

    if (onlyWeighMutants)
    {
        # If only weighting mutants, set mutations that are not mutants to a fixed
        # probability.

        mutationsTable <- mutationsTable %>%
            mutate(REAL_LENGTH_PROB_NORMAL = ifelse(MUTANT, REAL_LENGTH_PROB_NORMAL, 0.1),
                   REAL_LENGTH_PROB_MUTANT = ifelse(MUTANT, REAL_LENGTH_PROB_MUTANT, 0.1))
    }

    likelihoodRatio <-
        calc_likelihood_ratio_with_RL(M = mutationsTable$MUTANT,
                                      R = mutationsTable$DP,
                                      AF = mutationsTable$TUMOUR_AF,
                                      e = mutationsTable$BACKGROUND_AF,
                                      RL = mutationsTable$SIZE,
                                      RL_PROB_0 = mutationsTable$REAL_LENGTH_PROB_NORMAL,
                                      RL_PROB_1 = mutationsTable$REAL_LENGTH_PROB_MUTANT)

    likelihoodRatioListToTibble(likelihoodRatio)
}

calculateLikelihoodRatioForSampleWithoutSize <- function(mutationsTable, sizeTable)
{
    # Use depth = 1 for these calculations.

    mutationsTable <- mutationsTable %>%
        mutate(DP = 1)

    # Since not using sizes, this function is much smaller than the counterpart.

    likelihoodRatio <-
        calc_likelihood_ratio(M = mutationsTable$MUTANT,
                              R = mutationsTable$DP,
                              AF = mutationsTable$TUMOUR_AF,
                              e = mutationsTable$BACKGROUND_AF)

    likelihoodRatioListToTibble(likelihoodRatio)
}


##
# A single iteration of the body of doMain, as a function to assist with lapply.
#

singleIteration <- function(iteration, mutationsTable, sizeTable,
                            minFragmentLength, maxFragmentLength,
                            smooth, onlyWeighMutants,
                            .samplingSeed = NA)
{
    sampled <- mutationsTable
    if (iteration > 1)
    {
        if (!is.na(.samplingSeed))
        {
            set.seed(.samplingSeed)
        }

        sampled <- mutationsTable %>%
            slice_sample(n = nrow(mutationsTable), replace = TRUE)
    }

    usingSize <-
        calculateLikelihoodRatioForSampleWithSize(sampled, sizeTable,
                                                  minFragmentLength = minFragmentLength,
                                                  maxFragmentLength = maxFragmentLength,
                                                  smooth = smooth,
                                                  onlyWeighMutants = onlyWeighMutants) %>%
        mutate(USING_SIZE = TRUE)

    noSize <-
        calculateLikelihoodRatioForSampleWithoutSize(sampled, sizeTable) %>%
        mutate(USING_SIZE = FALSE)

    bind_rows(usingSize, noSize) %>%
        rename(INVAR_SCORE = LR, AF_P = P_MLE) %>%
        mutate(ITERATION = iteration)
}

##
# Function to return an INVAR scores table with headers but no rows.
# This can be the case where either there are no mutations
# after whole table filtering on entry to this function or when,
# for a subset of conditions, there are no rows left.
# This can optionally be an empty table for pre join (no SAMPLE_ID,
# PATIENT or PATIENT_MUTATION_BELONGS_TO) or post join (with these columns).
#

emptyInvarTable <- function(allColumns)
{
    assert_that(is.logical(allColumns), msg = "allColumns must be a logical")

    table <-
        tibble(ITERATION = integer(), USING_SIZE = logical(),
               LOCUS_NOISE.PASS = logical(), BOTH_STRANDS.PASS = logical(),
               OUTLIER.PASS = logical(), CONTAMINATION_RISK.PASS = logical(),
               INVAR_SCORE = double(), AF_P = double(),
               NULL_LIKELIHOOD = double(), ALTERNATIVE_LIKELIHOOD = double(),
               DP = integer(), MUTATION_SUM = integer(),
               IMAF = double(), SMOOTH = double(),
               OUTLIER_SUPPRESSION = double(), MUTANT_READS_PRESENT = logical())
     
    if (allColumns)
    {
        table <- table %>%
            add_column(SAMPLE_ID = character(), PATIENT = character(),
                       PATIENT_MUTATION_BELONGS_TO = character(),
                       .before = 1)
    }
    
    table
}

##
# Do the main processing for a single sample and patient mutation pair and
# outlier filter value.
#

doMain <- function(criteria, scriptArgs, mutationsTable, sizeTable, mc.set.seed = FALSE)
{
    assert_that(nrow(criteria) == 1, msg = "criteria expected to be exactly one row")

    # Calculate IMAFv2. Now have rows only for one sample and patient mutation
    # is for, so can use a slightly simpler function that returns a single value.

    IMAFv2 <- calculateIMAFv2(mutationsTable, bloodspot = scriptArgs$BLOODSPOT)

    # Apply standard filters to mutation table

    if (criteria$LOCUS_NOISE.PASS)
    {
        mutationsTable <- mutationsTable %>%
            filter(LOCUS_NOISE.PASS)
    }
    if (criteria$BOTH_STRANDS.PASS)
    {
        mutationsTable <- mutationsTable %>%
            filter(BOTH_STRANDS.PASS)
    }
    if (criteria$OUTLIER.PASS)
    {
        mutationsTable <- mutationsTable %>%
            filter(OUTLIER.PASS)
    }

    # Retune size characterisation table.

    sizeTable <- leaveOneOutFilter(mutationsTable, sizeTable)

    # tumour AF and size filters
    # for error classes where there are zero mutant reads, add one
    # just so that the background error estimate isn't zero

    mutationsTable <- mutationsTable %>%
        filter(TUMOUR_AF > 0 & SIZE > scriptArgs$MINIMUM_FRAGMENT_LENGTH & SIZE <= scriptArgs$MAXIMUM_FRAGMENT_LENGTH) %>%
        mutate(BACKGROUND_AF = ifelse(BACKGROUND_AF > 0, BACKGROUND_AF, 1 / BACKGROUND_DP))

    if (nrow(mutationsTable) == 0)
    {
        # If there is nothing left, return an empty INVAR scores table and print a warning.

        warning("No mutations left after filtering for LOCUS_NOISE.PASS = ", criteria$LOCUS_NOISE.PASS,
                ", BOTH_STRANDS.PASS = ", criteria$BOTH_STRANDS.PASS, ", OUTLIER.PASS = ", criteria$OUTLIER.PASS)

        return(emptyInvarTable(FALSE))
    }

    # determine whether there are any mutant reads in the table for calculation of ctDNA

    mutantReadsPresent <- any(mutationsTable$MUTANT)

    # The main method's filter on CONTAMINATION_RISK.PASS
    # implies that contaminationRisk is always only FALSE.

    contaminationRisk <- unique(mutationsTable$CONTAMINATION_RISK.PASS)
    assert_that(length(contaminationRisk) == 1, msg = "Have mix of CONTAMINATION_RISK.PASS flags in mutations table rows.")

    # There are two rows per molecule, so from the sorted table we can
    # slice out every other row to give a row per molecule.

    mutationsTable.perMolecule <- mutationsTable %>%
        arrange(UNIQUE_POS, SAMPLE_ID, SIZE) %>%
        slice(seq(1, n(), by = 2))

    # only count mutant reads once for accurate ctDNA quantification

    patientSpecific <- mutationsTable %>%
        distinct(PATIENT, PATIENT_MUTATION_BELONGS_TO) %>%
        summarise(PATIENT_SPECIFIC = PATIENT == PATIENT_MUTATION_BELONGS_TO)
    assert_that(nrow(patientSpecific) == 1, msg = "Have more than one patient + patient mutation belongs to pairing in file")

    iterations <- ifelse(patientSpecific$PATIENT_SPECIFIC, 1, 10)

    allIterations <-
        mclapply(1:iterations, singleIteration,
                 mutationsTable.perMolecule, sizeTable,
                 minFragmentLength = scriptArgs$MINIMUM_FRAGMENT_LENGTH,
                 maxFragmentLength = scriptArgs$MAXIMUM_FRAGMENT_LENGTH,
                 smooth = scriptArgs$SMOOTHING,
                 onlyWeighMutants = scriptArgs$ONLY_WEIGH_MUTANTS,
                 .samplingSeed = scriptArgs$SAMPLING_SEED,
                 mc.cores = scriptArgs$THREADS, mc.set.seed = mc.set.seed)

    mutationsTableSummary <- mutationsTable %>%
        summarise(DP = n(), MUTATION_SUM = sum(MUTANT))

    combinedResults <-
        bind_rows(allIterations) %>%
        full_join(criteria, by = character()) %>%
        full_join(mutationsTableSummary, by = character()) %>%
        mutate(IMAF = IMAFv2,
               SMOOTH = scriptArgs$SMOOTHING,
               OUTLIER_SUPPRESSION = scriptArgs$OUTLIER_SUPPRESSION,
               CONTAMINATION_RISK.PASS = contaminationRisk,
               MUTANT_READS_PRESENT = mutantReadsPresent)

    combinedResults
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    hasRNGSeed <- is.na(scriptArgs$SAMPLING_SEED)
    if (hasRNGSeed)
    {
        # See the help for mcparallel
        RNGkind("L'Ecuyer-CMRG")
        mc.reset.stream()
    }

    assert_that(file.exists(scriptArgs$MUTATIONS_TABLE_FILE),
                msg = str_c(scriptArgs$MUTATIONS_TABLE_FILE, " does not exist."))

    assert_that(file.exists(scriptArgs$SIZE_CHARACTERISATION_FILE),
                msg = str_c(scriptArgs$SIZE_CHARACTERISATION_FILE, " does not exist."))

    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    if (nrow(mutationsTable) == 0)
    {
        warning("Have no mutations in ", basename(scriptArgs$MUTATIONS_TABLE_FILE))
    }

    sizeTable <-
        readRDS(scriptArgs$SIZE_CHARACTERISATION_FILE) %>%
        filter(PATIENT_SPECIFIC) %>%
        select(MUTANT, SIZE, COUNT, TOTAL, PROPORTION)

    # For the time being, remove all rows that fail CONTAMINATION_RISK.PASS
    # Might want to revisit this decision, and possibly make it a condition
    # in the "crossing" table below.

    mutationsTable <- mutationsTable %>%
        filter(CONTAMINATION_RISK.PASS)

    if (nrow(mutationsTable) == 0)
    {
        warning("Have no mutations in ", basename(scriptArgs$MUTATIONS_TABLE_FILE), " after filtering for contamination.")

        # Create an empty table for saving.

        invarResultsTable <- emptyInvarTable(TRUE)
    }
    else
    {
        mutationsInfo <- mutationsTable %>%
            distinct(SAMPLE_ID, PATIENT, PATIENT_MUTATION_BELONGS_TO)

        # MutationsInfo should be 1x3
        assert_that(nrow(mutationsInfo) == 1, msg = "Do not have unique SAMPLE_ID, PATIENT, PATIENT_MUTATION_BELONGS_TO in mutations table file.")

        mutationsFileCheck <- mutationsInfo %>%
            mutate(MATCHING = SAMPLE_ID == scriptArgs$SAMPLE_ID & PATIENT_MUTATION_BELONGS_TO == scriptArgs$PATIENT_ID)

        assert_that(all(mutationsFileCheck$MATCHING),
                    msg = str_c("Mutations in", scriptArgs$MUTATIONS_TABLE_FILE, "do not belong to ",
                                scriptArgs$SAMPLE_ID, "and patient (mutation belongs to)",
                                scriptArgs$PATIENT_ID, sep = " "))

        # Create a table of all combinations of variable filter, then turn this into
        # a list of single row tibbles. lapply can then be used to work on every
        # sample + filter combination, and within that iterations can be run in parallel.

        allFilterCombinations <-
            crossing(OUTLIER.PASS = c(TRUE, FALSE),
                     LOCUS_NOISE.PASS = TRUE,
                     BOTH_STRANDS.PASS = TRUE)

        slicer <- function(n, table) { slice(table, n) }

        # False True True, or True True True
        allFilterCombinationList <-
            lapply(1:nrow(allFilterCombinations), slicer, allFilterCombinations)

        # Loop over list of files and allFilterCombinationList options
        invarResultsList <-
            lapply(allFilterCombinationList, doMain,
                   scriptArgs, mutationsTable, sizeTable,
                   mc.set.seed = hasRNGSeed)

        invarResultsTable <- bind_rows(invarResultsList)
        
        # If there a no rows at all, return an empty table.
        
        if (nrow(invarResultsTable) == 0)
        {
            warning("No rows in INVAR results table.")
            invarResultsTable <- emptyInvarTable(TRUE)
        }
        else
        {
            invarResultsTable <- invarResultsTable %>%
                full_join(mutationsInfo, by = character()) %>%
                select(SAMPLE_ID, PATIENT, PATIENT_MUTATION_BELONGS_TO,
                       ITERATION, USING_SIZE,
                       LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS, CONTAMINATION_RISK.PASS,
                       INVAR_SCORE, AF_P, NULL_LIKELIHOOD, ALTERNATIVE_LIKELIHOOD,
                       DP, MUTATION_SUM, IMAF, SMOOTH, OUTLIER_SUPPRESSION, MUTANT_READS_PRESENT) %>%
                arrange(SAMPLE_ID, PATIENT_MUTATION_BELONGS_TO,
                        ITERATION, USING_SIZE, LOCUS_NOISE.PASS, BOTH_STRANDS.PASS, OUTLIER.PASS)
        }
    }

    outputName <- str_c("invar_scores", makeSafeForFileName(scriptArgs$SAMPLE_ID),
                        makeSafeForFileName(scriptArgs$PATIENT_ID), "rds", sep = ".")

    saveRDS(invarResultsTable, outputName)
}

# Launch it.

invisible(main(parseOptions()))
