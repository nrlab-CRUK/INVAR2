##
# This script is essentially the "parse" function from TAPAS_functions.R
# plus the functions it (parse) uses.
# It's quite big enough as a script before going into all the processing
# that the function "TAPAS4" invokes after calling parse.
#
# Loads the mutations table and filters for on target mutations: those
# indicated in the patient file. Adds background error rates per mutation
# and labels each mutation as specific to a patient or not.
##

suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))


##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"

    options_list <- list(
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by createMutationsTable.R",
                    default=defaultMarker),
        make_option(c("--tumour-mutations"), type="character", metavar="file",
                    dest="TUMOUR_MUTATIONS_FILE", help="The source patient mutations file",
                    default=defaultMarker),
        make_option(c("--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("--error-rates"), type="character", metavar="file",
                    dest="ERROR_RATES_FILE", help="The error rates RDS file saved by offTargetErrorRate.R",
                    default=defaultMarker))

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
    testhome <- str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/')
    base <- str_c(testhome, 'createOnTargetMutationsTable/source/')

    list(
        MUTATIONS_TABLE_FILE = str_c(base, 'mutation_table.filtered.rds'),
        ERROR_RATES_FILE = str_c(base, 'error_rates.off_target.no_cosmic.rds'),
        TUMOUR_MUTATIONS_FILE = str_c(testhome, 'invar_source/PARADIGM_mutation_list_full_cohort_hg19.v2.csv'),
        LAYOUT_FILE = str_c(testhome, 'invar_source/combined.SLX_table_with_controls_031220.v2.csv')
    )
}


##
# From TAPAS_functions.R "combine_classes_in_rds"
#
# Complement the MUTATION_CLASS and TRINUCLEOTIDE columns for reference alleles
# 'A' and 'G'. 'T' and 'C' remain unchanged.
#
convertComplementaryMutations <- function(mutationTable)
{
    complement <- function(sequence)
    {
        chartr('ATCG', 'TAGC', sequence)
    }

    reverseComplement <- function(sequence)
    {
        stringi::stri_reverse(complement(sequence))
    }

    complementary <- function(base)
    {
        base == 'A' | base == 'G'
    }

    forward <- mutationTable %>%
        filter(!complementary(REF))

    reverse <- mutationTable %>%
        filter(complementary(REF))

    reverse <- reverse %>%
        mutate(TRINUCLEOTIDE = reverseComplement(TRINUCLEOTIDE),
               MUTATION_CLASS = complement(MUTATION_CLASS))

    bind_rows(forward, reverse) %>%
        arrange(SAMPLE_ID, CHROM, POS, REF, ALT, TRINUCLEOTIDE)
}


##
# From TAPAS_functions.R, a subpart of "parse"
#
# Classify mutations as being patient specific or not. Adds the TUMOUR_AF
# value from the tumour mutations table.
#
classifyForPatientSpecificity <- function(mutationTable, tumourMutationTable, layoutTable)
{
    # Patient specific is an inner join with the tumour mutation table by unique patient position.
    # Inner join combined the semi join with the addition of the TUMOUR_AF and MUTATION_CLASS
    # columns done in two steps in the original.

    tumourMutationTable.specific <- tumourMutationTable %>%
        select(UNIQUE_PATIENT_POS, TUMOUR_AF, MUTATION_CLASS, PATIENT)

    patientSpecific <- mutationTable %>%
        inner_join(tumourMutationTable.specific, by = 'UNIQUE_PATIENT_POS') %>%
        rename(PATIENT = PATIENT.x, PATIENT_MUTATION_BELONGS_TO = PATIENT.y)

    # Non-patient specific is basically the rows that do not match a patient specific record.
    # The TUMOUR_AF and MUTATION_CLASS values come from a unique position in the tumour
    # mutations table.

    tumourMutationTable.nonspecific <- tumourMutationTable %>%
        filter(!duplicated(UNIQUE_POS)) %>%
        select(UNIQUE_POS, TUMOUR_AF, MUTATION_CLASS, PATIENT)

    nonPatientSpecific <- mutationTable %>%
        anti_join(tumourMutationTable.specific, by = 'UNIQUE_PATIENT_POS') %>%
        inner_join(tumourMutationTable.nonspecific, by = 'UNIQUE_POS') %>%
        rename(PATIENT = PATIENT.x, PATIENT_MUTATION_BELONGS_TO = PATIENT.y)

    controlSamples <- layoutTable %>%
        filter(str_detect(CASE_OR_CONTROL, "control"))

    # ensure only loci that are interrogated in ptspec are looked for in the nonptspec

    patientSpecific <- patientSpecific %>%
        filter(UNIQUE_POS %in% nonPatientSpecific$UNIQUE_POS)

    nonPatientSpecific.cases <- nonPatientSpecific %>%
        filter(!SAMPLE_ID %in% controlSamples$SAMPLE_ID) %>%
        filter(UNIQUE_POS %in% patientSpecific$UNIQUE_POS)

    nonPatientSpecific.controls <- nonPatientSpecific %>%
        filter(SAMPLE_ID %in% controlSamples$SAMPLE_ID) %>%
        filter(UNIQUE_POS %in% patientSpecific$UNIQUE_POS)

    # Combine these three to provide a final mutation table with TUMOUR_AF values.

    bind_rows(patientSpecific, nonPatientSpecific.cases, nonPatientSpecific.controls)
}


##
# From TAPAS_functions.R "calculate.background_error"
#
calculateBackgroundError <- function(errorRatesList, layoutTable)
{
    # If the locus noise dataframe is completely clean - this means you have either a) set the locus noise threshold too low or b) not run enough control samples
    if (sum(errorRatesList$LOCUS_NOISE$MUTATION_SUM) == 0)
    {
        stop(str_c("LOCUS NOISE mutant sum is zero.",
                   "This is likely due to insufficient controls being run.",
                   "Either run more controls or override this with a higher locus noise filter (beware higher background rates).",
                   sep=' '))
    }

    thinLayoutTable <- layoutTable %>%
        select(SAMPLE_ID, CASE_OR_CONTROL)

    allErrorRates <- bind_rows(errorRatesList) %>%
        left_join(thinLayoutTable, by = 'SAMPLE_ID')

    backgroundError <- allErrorRates %>%
        group_by(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
        summarize(MUTATION_SUM = sum(MUTATION_SUM),
                  DP_SUM = sum(DP_SUM),
                  BACKGROUND_AF = MUTATION_SUM / DP_SUM,
                  .groups="drop")

    trinucleotideDepth <- allErrorRates %>%
        group_by(TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
        summarize(TRINUCLEOTIDE_DEPTH = sum(DP_SUM), .groups="drop")

    backgroundError2 <- backgroundError %>%
        filter(ALT != '.') %>%
        group_by(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
        inner_join(trinucleotideDepth, by = c('TRINUCLEOTIDE', 'CASE_OR_CONTROL', 'ERROR_RATE_TYPE')) %>%
        summarize(MUTATION_SUM_TOTAL = sum(MUTATION_SUM),
                  BACKGROUND_AF = MUTATION_SUM_TOTAL / TRINUCLEOTIDE_DEPTH,
                  TRINUCLEOTIDE_DEPTH = TRINUCLEOTIDE_DEPTH,
                  .groups="drop") %>%
        arrange(REF, ALT, CASE_OR_CONTROL, ERROR_RATE_TYPE)

    missingClasses <- missingErrorClasses(backgroundError2, trinucleotideDepth)

    bind_rows(backgroundError2, missingClasses) %>%
        arrange(REF, ALT, CASE_OR_CONTROL, ERROR_RATE_TYPE)
}


##
# From TAPAS_functions.R "add.missing_error_classes"
#
missingErrorClasses <- function(backgroundErrorTable, trinucleotideDepthTable)
{
    threeColumnJoin <- c('TRINUCLEOTIDE', 'CASE_OR_CONTROL', 'ERROR_RATE_TYPE')
    fiveColumnJoin <- c('REF', 'ALT', threeColumnJoin)

    uniqueClasses <- backgroundErrorTable %>%
        distinct(REF, ALT, TRINUCLEOTIDE)

    types = backgroundErrorTable %>%
        group_by(CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
        summarise(.groups = 'drop')

    allClasses <- uniqueClasses %>%
        crossing(types)

    presentClasses <- backgroundErrorTable %>%
        distinct(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE)

    missingClasses <- allClasses %>%
        anti_join(presentClasses, by = fiveColumnJoin) %>%
        mutate(MUTATION_SUM_TOTAL = 0, BACKGROUND_AF = 0.0) %>%
        left_join(trinucleotideDepthTable, by = threeColumnJoin)

    missingClasses
}

##
# From TAPAS_functions.R "combine_classes"
#
# Complement the REF, ALT and TRINUCLEOTIDE columns for reference alleles
# 'A' and 'G'. 'T' and 'C' remain unchanged.
#
convertComplementaryBackgroundError <- function(backgroundErrorTable)
{
    complement <- function(sequence)
    {
        chartr('ATCG', 'TAGC', sequence)
    }

    reverseComplement <- function(sequence)
    {
        stringi::stri_reverse(complement(sequence))
    }

    complementary <- function(base)
    {
        base == 'A' | base == 'G'
    }

    forward <- backgroundErrorTable %>%
        filter(!complementary(REF))

    reverse <- backgroundErrorTable %>%
        filter(complementary(REF))

    reverse <- reverse %>%
        mutate(TRINUCLEOTIDE = reverseComplement(TRINUCLEOTIDE),
               REF = complement(REF),
               ALT = complement(ALT))

    bind_rows(forward, reverse) %>%
        arrange(REF, ALT, TRINUCLEOTIDE)
}



##
# From TAPAS_functions.R "combine_classes"
#
# "combine complementary classes and trinucleotide contexts i.e. A>T etc."
# In effect, reverse the strand and change 'A' to 'T' and 'G' to 'C'.
# C & T are not changed.
#
errorTableComplementaryClasses <- function(backgroundErrorTable)
{
    complement <- backgroundErrorTable %>%
        convertComplementaryBackgroundError() %>%
        group_by(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
        summarise(MUTATION_SUM_TOTAL = sum(MUTATION_SUM_TOTAL),
                  TRINUCLEOTIDE_DEPTH = sum(TRINUCLEOTIDE_DEPTH),
                  BACKGROUND_AF = MUTATION_SUM_TOTAL / TRINUCLEOTIDE_DEPTH,
                  .groups = "drop")

    complement
}

##
# Referred functions are from TAPAS_functions.R
#

addPatientAndBackgroundColumns <- function(mutationTable, tumourMutationTable, layoutTable, backgroundErrorTable)
{
    # This bit from "annotate_with_SLX_table", plus the additional column

    mutationTable <- mutationTable %>%
        left_join(layoutTable, by = 'SAMPLE_ID') %>%
        mutate(UNIQUE_PATIENT_POS = str_c(PATIENT, UNIQUE_POS, sep='_'))

    # This from "ignore_unexpected_variants". After this point, we want to
    # remove MUTATION and UNIQUE_ALT from the mutation table as their values
    # will come from joining with the tumour mutation table
    # (in classifyForPatientSpecificity).

    mutationTable <- mutationTable %>%
        filter(AF == 0 | UNIQUE_ALT %in% tumourMutationTable$UNIQUE_ALT)

    # Remainder from "parse"

    # Add patient specific indicator and the TUMOUR_AF value from the tumour
    # mutations table.
    # Also convert to complementary strand for A and G reference alleles.

    mutationTable.withPatient <- mutationTable %>%
        classifyForPatientSpecificity(tumourMutationTable, layoutTable) %>%
        convertComplementaryMutations()

    # Add background error rates.

    backgroundErrorTable.forJoin <- backgroundErrorTable %>%
        mutate(MUTATION_CLASS = str_c(REF, ALT, sep = '/')) %>%
        rename(BACKGROUND_MUTATION_SUM = MUTATION_SUM_TOTAL,
               BACKGROUND_DP = TRINUCLEOTIDE_DEPTH) %>%
        filter(CASE_OR_CONTROL == 'case' & ERROR_RATE_TYPE == 'locus_noise.both_strands') %>%
        select(TRINUCLEOTIDE, MUTATION_CLASS, starts_with('BACKGROUND_'))

    # Add background columns to mutation table.

    mutationTable.withPatientAndBackground <- mutationTable.withPatient %>%
        inner_join(backgroundErrorTable.forJoin, by = c('TRINUCLEOTIDE', 'MUTATION_CLASS'))

    mutationTable.withPatientAndBackground
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    assert_that(file.exists(scriptArgs$TUMOUR_MUTATIONS_FILE), msg = str_c(scriptArgs$TUMOUR_MUTATIONS_FILE, " does not exist."))
    assert_that(file.exists(scriptArgs$LAYOUT_FILE), msg = str_c(scriptArgs$LAYOUT_FILE, " does not exist."))
    assert_that(file.exists(scriptArgs$MUTATIONS_TABLE_FILE), msg = str_c(scriptArgs$MUTATIONS_TABLE_FILE, " does not exist."))
    assert_that(file.exists(scriptArgs$ERROR_RATES_FILE), msg = str_c(scriptArgs$ERROR_RATES_FILE, " does not exist."))

    tumourMutationTable <-
        loadTumourMutationsTable(scriptArgs$TUMOUR_MUTATIONS_FILE) %>%
        select(PATIENT, REF, ALT, TUMOUR_AF, MUTATION_CLASS, UNIQUE_POS, UNIQUE_ALT, UNIQUE_PATIENT_POS)

    layoutTable <-
        loadLayoutTable(scriptArgs$LAYOUT_FILE) %>%
        select(SAMPLE_ID, PATIENT, CASE_OR_CONTROL)

    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE)

    if (nrow(mutationTable) == 0)
    {
        stop("There are no mutations in the mutations table (", scriptArgs$MUTATIONS_TABLE_FILE, ").")
    }

    mutationTable <- mutationTable %>%
        filter(ON_TARGET) %>%
        addMutationTableDerivedColumns()

    if (nrow(mutationTable) == 0)
    {
        stop("There are no on-target mutations in the mutations table.")
    }

    errorRatesList <- readRDS(scriptArgs$ERROR_RATES_FILE)

    for (v in names(errorRatesList))
    {
        errorRatesList[[v]] = errorRatesList[[v]] %>%
            mutate(ERROR_RATE_TYPE = str_to_lower(v))
    }

    backgroundErrorTable <- errorRatesList %>%
        calculateBackgroundError(layoutTable) %>%
        errorTableComplementaryClasses()

    mutationTable.withPatientAndBackground <-
        addPatientAndBackgroundColumns(mutationTable, tumourMutationTable, layoutTable, backgroundErrorTable)

    if (nrow(mutationTable.withPatientAndBackground) == 0)
    {
        stop("There are no on-target mutations after filtering for background.")
    }

    backgroundErrorTable %>%
        saveRDS('background_error_rates.rds')

    mutationTable.withPatientAndBackground %>%
        removeMutationTableDerivedColumns() %>%
        arrangeMutationTableForExport() %>%
        saveRDS('mutation_table.on_target.all.rds')
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
