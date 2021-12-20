# suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))


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
                    default=defaultMarker),
        make_option(c("--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."),
        make_option(c("--control-proportion"), type="double", metavar="num",
                    dest="CONTROL_PROPORTION", help="Blacklist loci that have signal in >30% of the nonptspec samples",
                    default=0.1),
        make_option(c("--max-background-af"), type="double", metavar="num",
                    dest="MAX_BACKGROUND_AF", help="Filter loci with a background AF in controls greater than this value",
                    default=0.01))

    opts <- OptionParser(option_list=options_list, usage="%prog [options] <mutation file>") %>%
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
        MUTATIONS_TABLE_FILE = 'mutations/mutation_table.filtered.rds',
        TUMOUR_MUTATIONS_FILE = 'source_files/PARADIGM_mutation_list_full_cohort_hg19.csv',
        LAYOUT_FILE = 'source_files/combined.SLX_table_with_controls_031220.csv',
        ERROR_RATES_FILE = 'off_target/mutation_table.error_rates.no_cosmic.rds',
        CONTROL_PROPORTION = 0.1,
        MAX_BACKGROUND_AF = 0.01,
        BLOODSPOT = FALSE
    )
}


##
# Loading functions.
#


# Load the patient specific tumour mutations file

loadTumourMutationsTable <- function(tumourMutationsFile)
{
    read_csv(tumourMutationsFile, col_types = 'ciccd', show_col_types = FALSE) %>%
        select(-contains('uniq')) %>%
        rename_with(str_to_upper) %>%
        rename(CHROM = CHR) %>%
        mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':'),
               MUTATION = str_c(REF, ALT, sep='/'),
               UNIQUE_ALT = str_c(UNIQUE_POS, MUTATION, sep='_'),
               UNIQUE_PATIENT_POS = str_c(PATIENT, UNIQUE_POS, sep='_')) %>%
        select(PATIENT, REF, ALT, TUMOUR_AF, UNIQUE_POS, UNIQUE_ALT, UNIQUE_PATIENT_POS)
}

# Read the layout file and extract unique pool id and barcode pairs.
# Also only keep the columns this script needs to use.

loadLayoutFile <- function(layoutFile)
{
    suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE)) %>%
        rename_with(str_to_upper) %>%
        mutate(POOL_BARCODE = str_c(SLX_ID, str_replace(BARCODE, '-', '_'), sep = '_')) %>%
        select(STUDY, SAMPLE_NAME, PATIENT, SAMPLE_TYPE, CASE_OR_CONTROL, INPUT_INTO_LIBRARY_NG, POOL_BARCODE)
}

addDerivedColumns <- function(mutationTable)
{
    mutationTable %>%
        mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':'),
               POOL_BARCODE = str_c(POOL, BARCODE, sep='_'),
               MUTATION = str_c(REF, ALT, sep='/'),
               UNIQUE_ALT = str_c(UNIQUE_POS, MUTATION, sep='_'))
}

##
# From TAPAS_functions.R, originally "annotate_with_locus_error_rate"
#
# Originally just annotate, this method creates the errorRateTable and returns.
#
# Locus error rate = overall background error rate per locus, aggregated across control samples
createErrorRateTable <- function(mutationTable,
                                 layoutTable,
                                 proportion_of_controls,
                                 max_background_mean_AF,
                                 is.blood_spot)
{
    if (is.blood_spot)
    {
        message("sWGS/blood spot mutationTable, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
        max_background_mean_AF <- 1
    }

    errorRateTable <- mutationTable %>%
        filter(POOL_BARCODE %in% layoutTable.cases$POOL_BARCODE) %>%
        mutate(HAS_SIGNAL = ifelse(ALT_F + ALT_R > 0, POOL_BARCODE, NA)) %>%
        group_by(UNIQUE_POS, TRINUCLEOTIDE) %>%
        summarize(MUT_SUM = sum(ALT_F) + sum(ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF = MUT_SUM / DP_SUM,
                  N_SAMPLES = n_distinct(POOL_BARCODE),
                  N_SAMPLES_WITH_SIGNAL = n_distinct(HAS_SIGNAL, na.rm = TRUE),
                  .groups = 'drop')

    errorRateTable %>%
        mutate(PROPORTION = N_SAMPLES_WITH_SIGNAL / N_SAMPLES,
               LOCUS_NOISE.PASS = PROPORTION < proportion_of_controls & BACKGROUND_AF < max_background_mean_AF)
}

##
# From TAPAS_functions.R "combine_classes_in_rds"
#
# Complement the MUTATION and TRINUCLEOTIDE columns for reference alleles
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

    mutationTable %>%
        mutate(TRINUCLEOTIDE = ifelse(complementary(REF), reverseComplement(TRINUCLEOTIDE), TRINUCLEOTIDE),
               MUTATION = ifelse(complementary(REF), complement(MUTATION), MUTATION))
}


##
# From TAPAS_functions.R, a subpart of "parse"
#
# Classify mutations as being patient specific or not. Adds the TUMOUR_AF
# value from the tumour mutations table.
#
classifyForPatientSpecificity <- function(mutationTable, tumourMutationTable, layoutTable,
                                          skipPatientFromBackground = NULL)
{
    # Patient specific is an inner join with the tumour mutation table by unique patient position.
    # Inner join combined the semi join with the addition of the TUMOUR_AF column done in two
    # steps in the original.

    tumourMutationTable.specific <- tumourMutationTable %>%
        select(UNIQUE_PATIENT_POS, TUMOUR_AF)

    patientSpecific <- mutationTable %>%
        inner_join(tumourMutationTable.specific, by = 'UNIQUE_PATIENT_POS')

    if (!is.null(skipPatientFromBackground))
    {
        patientSpecific <- patientSpecific %>%
            filter(!str_detect(skip_patient_from_background))
    }

    # Non-patient specific is basically the rows that do not match a patient specific record.
    # The TUMOUR_AF value comes from a unique position in the tumour mutations table.

    tumourMutationTable.nonspecific <- tumourMutationTable %>%
        filter(!duplicated(UNIQUE_POS)) %>%
        select(UNIQUE_POS, TUMOUR_AF)

    nonPatientSpecific <- mutationTable %>%
        anti_join(tumourMutationTable.specific, by = 'UNIQUE_PATIENT_POS') %>%
        inner_join(tumourMutationTable.nonspecific, by = 'UNIQUE_POS')

    controlSamples <- layoutTable %>%
        filter(str_detect(CASE_OR_CONTROL, "control"))

    # ensure only loci that are interrogated in ptspec are looked for in the nonptspec

    patientSpecific <- patientSpecific %>%
        filter(UNIQUE_POS %in% nonPatientSpecific$UNIQUE_POS) %>%
        mutate(PATIENT_SPECIFIC = TRUE)

    nonPatientSpecific.cases <- nonPatientSpecific %>%
        filter(!POOL_BARCODE %in% controlSamples$POOL_BARCODE) %>%
        filter(UNIQUE_POS %in% patientSpecific$UNIQUE_POS) %>%
        mutate(PATIENT_SPECIFIC = FALSE)

    nonPatientSpecific.controls <- nonPatientSpecific %>%
        filter(POOL_BARCODE %in% controlSamples$POOL_BARCODE) %>%
        filter(UNIQUE_POS %in% patientSpecific$UNIQUE_POS) %>%
        mutate(PATIENT_SPECIFIC = FALSE)

    # Combine these three to provide a final mutation table with TUMOUR_AF values.

    bind_rows(patientSpecific, nonPatientSpecific.cases, nonPatientSpecific.controls)
}


##
# From TAPAS_functions.R "calculate.background_error"
# PPC = control sample naming in NR lab
#
calculateBackgroundError <- function(errorRatesList, layoutTable, exclude_PPC = FALSE)
{
    # If the locus noise dataframe is completely clean - this means you have either a) set the locus noise threshold too low or b) not run enough control samples
    if (sum(errorRatesList$LOCUS_NOISE_PASS$MUT_SUM) == 0)
    {
        stop(str_c("LOCUS NOISE mutant sum is zero.",
                   "This is likely due to insufficient controls being run.",
                   "Either run more controls or override this with a higher locus noise filter (beware higher background rates).",
                   sep='\n'))
    }

    thinLayoutTable <- layoutTable %>%
        select(SAMPLE_NAME, CASE_OR_CONTROL, POOL_BARCODE)

    allErrorRates <- bind_rows(errorRatesList) %>%
        mutate(POOL_BARCODE = str_c(POOL, BARCODE, sep='_')) %>%
        left_join(thinLayoutTable, by = 'POOL_BARCODE')

    if (exclude_PPC)
    {
        allErrorRates <- allErrorRates %>%
            filter(!str_detect(SAMPLE_NAME, "PPC"))
    }

    backgroundError <- allErrorRates %>%
        group_by(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
        summarize(MUT_SUM = sum(MUT_SUM),
                  DP_SUM = sum(DP_SUM),
                  BACKGROUND_AF = MUT_SUM / DP_SUM,
                  .groups="drop")

    trinucleotideDepth <- backgroundError %>%
        group_by(TRINUCLEOTIDE) %>%
        summarize(DP = sum(DP_SUM), .groups="drop")

    totalDP.byType <- backgroundError %>%
        group_by(TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
        summarize(GROUP_DP = sum(DP_SUM), .groups="drop")

    backgroundError2 <- backgroundError %>%
        filter(ALT != '.') %>%
        group_by(REF, ALT, TRINUCLEOTIDE, CASE_OR_CONTROL, ERROR_RATE_TYPE) %>%
        inner_join(totalDP.byType, by = c('TRINUCLEOTIDE', 'CASE_OR_CONTROL', 'ERROR_RATE_TYPE'))
        summarize(GROUP_DP = GROUP_DP,
                  MUT_SUM_TOTAL = sum(MUT_SUM),
                  BACKGROUND_AF = MUT_SUM_TOTAL / GROUP_DP,
                  .groups="drop")

    backgroundError2
}


##
# Referred functions are from TAPAS_functions.R
#

processTables <- function(mutationTable, tumourMutationTable, layoutTable, errorRatesList,
                          skipPatientFromBackground = NULL)
{
    # This bit from "annotate_with_SLX_table", plus the additional column

    mutationTable <- mutationTable %>%
        left_join(layoutTable, by = 'POOL_BARCODE') %>%
        mutate(UNIQUE_PATIENT_POS = str_c(PATIENT, UNIQUE_POS, sep='_'))

    # This from "ignore_unexpected_variants"

    mutationTable <- mutationTable %>%
        filter(AF == 0 | UNIQUE_ALT %in% tumourMutationTable$UNIQUE_ALT)

    # Remainder from "parse"

    # Add patient specific indicator and the TUMOUR_AF value from the tumour
    # mutations table.
    # Also convert to complementary strand for A and G reference alleles.

    mutationTable.withPatient <- mutationTable %>%
        classifyForPatientSpecificity(tumourMutationTable, layoutTable, skipPatientFromBackground) %>%
        convertComplementaryMutations()

    # Calculate error rates.

    # backgroundErrorRates <-calculateBackgroundError(errorRatesList, layoutTable)

    mutationTable.withPatient
}


removeDerivedColums <- function(mutationTable)
{
    mutationTable %>%
        select(-any_of(c('UNIQUE_POS', 'POOL_BARCODE', 'UNIQUE_ID')))
}

saveRDSandTSV <- function(t, file)
{
    saveRDS(t, file)

    tsv <- str_replace(file, "\\.rds$", ".tsv")

    if (tsv != file)
    {
        write_tsv(t, tsv)
    }

    t
}

cleanAndSaveRDSandTSV <- function(t, file)
{
    t %>%
        removeDerivedColums() %>%
        saveRDSandTSV(file)

    t
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    tumourMutationTable <- loadTumourMutationsTable(scriptArgs$TUMOUR_MUTATIONS_FILE)

    layoutTable <- loadLayoutFile(scriptArgs$LAYOUT_FILE)

    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        filter(ON_TARGET) %>%
        addDerivedColumns()

    errorRatesList <- readRDS(scriptArgs$ERROR_RATES_FILE)

    for (v in names(errorRatesList))
    {
        errorRatesList[[v]] = errorRatesList[[v]] %>%
            mutate(ERROR_RATE_TYPE = str_to_lower(v))
    }

    mutationTable.combined <- processTables(mutationTable, tumourMutationTable, layoutTable, errorRatesList)

    mutationTable.combined
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789') {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
