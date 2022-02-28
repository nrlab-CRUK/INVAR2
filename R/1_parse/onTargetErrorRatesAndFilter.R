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
        make_option(c("--study"), type="character", metavar="string",
                    dest="STUDY", help="The study name",
                    default=defaultMarker),
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by createOnTargetMutationsTable.R",
                    default=defaultMarker),
        make_option(c("--cosmic-threshold"), type="integer", metavar="int",
                    dest="COSMIC_THRESHOLD", help="Loci with >0 entries in COSMIC are considered as COSMIC mutations",
                    default=0L),
        make_option(c("--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."),
        make_option(c("--control-proportion"), type="double", metavar="num",
                    dest="CONTROL_PROPORTION", help="Blacklist loci that have signal in >30% of the nonptspec samples",
                    default=0.1),
        make_option(c("--max-background-allele-frequency"), type="double", metavar="num",
                    dest="MAX_BACKGROUND_ALLELE_FREQUENCY", help="Filter loci with a background allele frequency in controls greater than this value",
                    default=0.01),
        make_option(c("--allele-frequency-threshold"), type="double", metavar="num",
                    dest="ALLELE_FREQUENCY_THRESHOLD", help="Maximum allele frequency value for acceptable samples",
                    default=0.01))

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
    base <- str_c(testhome, 'onTargetErrorRatesAndFilter/source/')

    list(
        MUTATIONS_TABLE_FILE = str_c(base, 'mutation_table.on_target.all.rds'),
        STUDY = 'PARADIGM',
        CONTROL_PROPORTION = 0.1,
        MAX_BACKGROUND_ALLELE_FREQUENCY = 0.01,
        ALLELE_FREQUENCY_THRESHOLD = 0.01,
        BLOODSPOT = FALSE,
        COSMIC_THRESHOLD = 0L
    )
}


##
# From TAPAS_functions.R, originally "annotate_with_locus_error_rate"
#
# Originally just annotate, this method creates the lociErrorRateTable and returns.
#
# Locus error rate = overall background error rate per locus, aggregated across control samples
#
createLociErrorRateTable <- function(mutationTable,
                                     proportion_of_controls,
                                     max_background_mean_AF,
                                     is.blood_spot)
{
    assert_that(is.number(proportion_of_controls), msg = "proportion_of_controls must be a number.")
    assert_that(is.number(max_background_mean_AF), msg = "max_background_mean_AF must be a number.")
    assert_that(is.flag(is.blood_spot), msg = "is.blood_spot must be a logical.")

    if (is.blood_spot)
    {
        message("sWGS/blood spot mutationTable, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
        max_background_mean_AF <- 1
    }

    errorRateTable <- mutationTable %>%
        filter(!PATIENT_SPECIFIC & CASE_OR_CONTROL == 'case') %>%
        mutate(HAS_SIGNAL = ifelse(ALT_F + ALT_R > 0, SAMPLE_ID, NA)) %>%
        group_by(UNIQUE_POS, MUTATION_CLASS, TRINUCLEOTIDE, PATIENT_MUTATION_BELONGS_TO, COSMIC) %>%
        summarize(MUTATION_SUM = sum(ALT_F + ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF = MUTATION_SUM / DP_SUM,
                  N_SAMPLES = n_distinct(SAMPLE_ID),
                  N_SAMPLES_WITH_SIGNAL = n_distinct(HAS_SIGNAL, na.rm = TRUE),
                  .groups = 'drop') %>%
        separate(UNIQUE_POS, sep = ':', into = c('CHROM', 'POS'), remove = FALSE) %>%
        mutate(POS = as.integer(POS)) %>%
        select(UNIQUE_POS, CHROM, POS, MUTATION_CLASS, TRINUCLEOTIDE, PATIENT_MUTATION_BELONGS_TO, COSMIC,
               BACKGROUND_AF, MUTATION_SUM, DP_SUM, N_SAMPLES, N_SAMPLES_WITH_SIGNAL) %>%
        mutate(LOCUS_NOISE.PASS = (N_SAMPLES_WITH_SIGNAL / N_SAMPLES) < proportion_of_controls &
                                  BACKGROUND_AF < max_background_mean_AF)

    errorRateTable
}


# identify samples with ctDNA >1% so that they are flagged
# later in the pipeline they are not used as patient-controls because of a
# potential risk of contamination between pt-spec and non-pt-spec loci
getContaminatedSamples <- function(mutationTable, afThreshold)
{
    assert_that(is.number(afThreshold), msg = "afThreshold must be a number.")

    alleleFrequencyTable <- mutationTable %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS.PASS) %>%
        group_by(SAMPLE_ID, PATIENT_MUTATION_BELONGS_TO, PATIENT_SPECIFIC) %>%
        summarise(MUTATION_SUM = sum(ALT_F + ALT_R),
                  DP = sum(DP),
                  AF = MUTATION_SUM / DP,
                  .groups = 'drop') %>%
        rename(PATIENT = PATIENT_MUTATION_BELONGS_TO)

    ## any correlation between sample and it's nonptspec and ptspec AFs?

    patientSpecific = alleleFrequencyTable %>%
        filter(PATIENT_SPECIFIC)

    nonPatientSpecific = alleleFrequencyTable %>%
        filter(!PATIENT_SPECIFIC)

    test <- inner_join(patientSpecific, nonPatientSpecific, by = c('SAMPLE_ID'))

    lowSamples <- test %>%
        filter(AF.x < afThreshold)

    print(cor.test(lowSamples$AF.x, lowSamples$AF.y))

    doNotUse <- alleleFrequencyTable %>%
        filter(PATIENT_SPECIFIC & AF > afThreshold)

    doNotUse %>%
        distinct(SAMPLE_ID)
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    lociErrorRateTable <- mutationTable %>%
        createLociErrorRateTable(proportion_of_controls = scriptArgs$CONTROL_PROPORTION,
                                 max_background_mean_AF = scriptArgs$MAX_BACKGROUND_ALLELE_FREQUENCY,
                                 is.blood_spot = scriptArgs$BLOODSPOT)

    lociErrorRateTable %>%
        select(-UNIQUE_POS) %>%
        saveRDS('locus_error_rates.on_target.rds')

    lociErrorRateNoisePass <- lociErrorRateTable %>%
        filter(LOCUS_NOISE.PASS)

    mutationTable.filtered <- mutationTable %>%
        filter(PATIENT_SPECIFIC | COSMIC_MUTATIONS <= scriptArgs$COSMIC_THRESHOLD) %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% lociErrorRateNoisePass$UNIQUE_POS,
               BOTH_STRANDS.PASS = ALT_F > 0 & ALT_R > 0 | AF == 0)

    contaminatedSamples <- getContaminatedSamples(mutationTable.filtered, scriptArgs$ALLELE_FREQUENCY_THRESHOLD)

    mutationTable.filtered <- mutationTable.filtered %>%
        mutate(CONTAMINATION_RISK.PASS = !SAMPLE_ID %in% contaminatedSamples$SAMPLE_ID)

    mutationTable.filtered %>%
        removeMutationTableDerivedColumns() %>%
        arrangeMutationTableForExport() %>%
        saveRDS('mutation_table.on_target.rds')
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
