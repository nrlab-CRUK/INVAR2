suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))


##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"

    options_list <- list(
        make_option(c("-t", "--tapas"), type="character", metavar="string",
                    dest="TAPAS_SETTING", help="The TAPAS setting string",
                    default=defaultMarker),
        make_option(c("-m", "--mutations-list"), type="character", metavar="file",
                    dest="MUTATIONS_BED_FILE", help="The source mutations list file",
                    default=defaultMarker),
        make_option(c("-l", "--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("-b", "--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."))

    opts <- OptionParser(option_list=options_list, usage="%prog [options] <mutation file>") %>%
        parse_args(positional_arguments = TRUE)

    missingRequired <- which(opts$options == defaultMarker)
    if (any(missingRequired))
    {
        missing <- names(opts$options)[missingRequired]
        stop(str_c("ERROR: One or more required arguments have not been provided: ", str_c(missing, collapse = ' ')))
    }

    if (length(opts$args) != 1)
    {
        stop("The mutation file has not been provided.")
    }

    scriptOptions <- list(MUTATIONS_FILE = opts$args[1])
    for (v in names(opts$options))
    {
        scriptOptions[v] = opts$options[v]
    }

    scriptOptions
}

# Test options for my (Rich) local set up in RStudio.

richTestOptions <- function()
{
    #mutationsFile <- 'mutations/SLX-19722_SXTLI005.mutations.tsv'
    list(
        MUTATIONS_FILE = '/mnt/scratchb/bioinformatics/bowers01/invar_emma/EMMA/output_gz/PARADIGM.f0.9_s2.BQ_20.MQ_40.combined.final.ann.tsv',
        TAPAS_SETTING = 'f0.9_s2.BQ_30.MQ_60',
        MUTATIONS_BED_FILE = 'bed/PARADIGM_mutation_list_full_cohort_hg19.bed',
        LAYOUT_FILE = 'bed/combined.SLX_table_with_controls_031220.csv'
    )
}


##
# Taken from TAPAS_functions.R
# Multiallelic blacklisting.
#

# Filter out rows where the MSQB value is greater than the threshold.

blacklist.MQSB <- function(mutationTable, individual_MQSB_threshold)
{
    message("number of unique loci (pre MQSB): ", length(unique(mutationTable$UNIQUE_POS)))

    message("applying an individiual MQSB threshold of ", individual_MQSB_threshold)

    filteredMutationTable <- mutationTable %>%
        filter(MQSB > individual_MQSB_threshold)

    message("number of unique loci (post MQSB): ", length(unique(filteredMutationTable$UNIQUE_POS)))

    filteredMutationTable
}

# Create a list of loci that have more than the given number of alleles.

createMultiallelicBlacklist <- function(mutationTable, n_alt_alleles_threshold = 3, minor_alt_allele_threshold = 2)
{
    # Filter for rows with positive AF and
    # determine the number of alt alleles per UNIQUE_POS
    multiallelic <- mutationTable %>%
        filter(AF > 0) %>%
        mutate(MUT_SUM = ALT_F + ALT_R) %>%
        select(UNIQUE_POS, ALT, MUT_SUM) %>%
        arrange(UNIQUE_POS)

    # Identify loci with >1 alt alleles (e.g. both A>C and A>T at a position)
    duplicatePositions <- multiallelic %>%
        filter(duplicated(UNIQUE_POS))

    multiallelic <- multiallelic %>%
        filter(UNIQUE_POS %in% duplicatePositions)

    # Create a list of loci that have more than one allele.

    multiallelic_blacklist <- tibble('UNIQUE_POS', .rows = 0)

    if (nrow(multiallelic) > 0)
    {
        multiallelic_blacklist <- multiallelic %>%
            group_by(UNIQUE_POS) %>%
            mutate(N_ALT_ALLELES = n(), MIN = min(MUT_SUM), MAX = max(MUT_SUM)) %>%
            filter(N_ALT_ALLELES >= n_alt_alleles_threshold &
                   MIN >= minor_alt_allele_threshold &
                   MAX >= minor_alt_allele_threshold) %>%
            select(UNIQUE_POS) %>%
            distinct(UNIQUE_POS)
    }

    multiallelic_blacklist
}

# Tidy multiallelic sites. Filter out those that have more than the requested
# number of alleles.

blacklist.multiallelic <- function(mutationTable,
                                   n_alt_alleles_threshold = 3, minor_alt_allele_threshold = 2,
                                   blacklistFile = NULL)
{
    message("starting rows before multiallelic blacklist: ", nrow(mutationTable))

    multiallelic_blacklist <- createMultiallelicBlacklist(mutationTable, n_alt_alleles_threshold, minor_alt_allele_threshold)

    message("length of multiallelic blacklist is ", nrow(multiallelic_blacklist))

    # save the list of multiallelic loci for reference
    if (!is.null(blacklistFile))
    {
        write_tsv(multiallelic_blacklist, blacklistFile, col_names = FALSE, quote = 'none')
    }

    filteredMutationTable <- mutationTable

    if (nrow(multiallelic_blacklist) > 0)
    {
        filteredMutationTable <- mutationTable %>%
            filter(!UNIQUE_POS %in% multiallelic_blacklist)
    }

    message("rows after multiallelic blacklist filtering: ", nrow(filteredMutationTable))

    filteredMutationTable
}


##
# Loading functions.
#


# Load the mutations table from file and filter.

load.mutations.table <- function(mutationsFile, patientSpecificFile, tapasSetting,
                                 cosmic_threshold = 0, max_DP = 1500, min_ref_DP = 5,
                                 individual_MQSB_threshold = 0.01,
                                 n_alt_alleles_threshold = 3, minor_alt_allele_threshold = 2)
{
    mutationTableRDSFile <- str_c(tapasSetting, '.rds')

    if (file.exists(mutationTableRDSFile))
    {
        # When run in a pipeline, this will never exist.
        # In development, it might and can save some time to use it.
        warning("Read mutations table from previously saved ", mutationTableRDSFile)
        return(readRDS(mutationTableRDSFile))
    }

    message("Reading in table from ", mutationsFile)

    patient_specific <-
        read_tsv(patientSpecificFile,
                 col_names = c('CHROM', 'START', 'END', 'REF', 'ALT'),
                 col_types = 'ciicc') %>%
        mutate(UNIQUE_POS = str_c(CHROM, END, sep=':'))

    mutationTable <- read_tsv(mutationsFile, col_types = 'ciccicdddddcciidc')

    # Remove soft-masked repeats, identified by lowercase.
    # Add columns of derived values and combined identifiers.
    mutationTable <- mutationTable %>%
        filter(!(str_detect(REF, '[acgt]') | str_detect(ALT, '[acgt]'))) %>%
        mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':')) %>%
        mutate(SLX_BARCODE = str_c(SLX, BARCODE, sep='_')) %>%
        mutate(UNIQUE_ID = str_c(UNIQUE_POS, SLX_BARCODE, sep='_')) %>%
        mutate(AF = (ALT_F + ALT_R) / DP) %>%
        mutate(COSMIC = COSMIC_MUTATIONS > cosmic_threshold) %>%
        mutate(SNP = `1KG_AF` > 0) %>%
        mutate(ON_TARGET = UNIQUE_POS %in% patient_specific$UNIQUE_POS)


    # save on.target pre-filtering
    mutationTable.prefilter <- mutationTable %>%
        filter(ON_TARGET & nchar(ALT) == 1 & nchar(REF) == 1)

    saveRDS(mutationTable.prefilter, file = str_c(tapasSetting, ".on_target.prefilter.rds"))

    # Apply filters to blacklist loci
    mutationTable <- mutationTable %>%
        filter(DP < max_DP & !SNP & REF_R + REF_F >= min_ref_DP & (ON_TARGET | !COSMIC))

    mutationTable <- blacklist.MQSB(mutationTable, individual_MQSB_threshold)
    mutationTable <- blacklist.multiallelic(mutationTable, n_alt_alleles_threshold, minor_alt_allele_threshold)

    saveRDS(mutationTable, mutationTableRDSFile)

    mutationTable
}

# Read the layout file and extract unique pool id and barcode pairs.

load.patient.samples <- function(layoutFile)
{
    read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE) %>%
    filter(case_or_control == "case") %>%
    mutate(SLX_BARCODE = str_c(SLX_ID, str_replace(barcode, '-', '_'), sep = '_')) %>%
    select(distinct(SLX_BARCODE))
}

# Filter rows that are on target.

filter.for.ontarget <- function(mutationTable)
{
    mutationTable %>% filter(ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)
}

##
# From TAPAS_functions.R
#
#' Annotate with locus error rate
#' Locus error rate = overall background error rate per locus, aggregated across control samples
annotate_with_locus_error_rate <- function(mutationTable,
                                           patientSamples,
                                           proportion_of_controls = 0.1,
                                           max_background_mean_AF = 0.01,
                                           is.blood_spot = FALSE,
                                           on_target = TRUE,
                                           errorRateFile = NULL)
{

    if (is.blood_spot)
    {
        message("sWGS/blood spot mutationTable, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
        max_background_mean_AF <- 1
    }

    errorRateTable <- NULL

    if (!is.null(errorRateFile) & file.exists(errorRateFile))
    {
        # To speed up testing.
        warning("Reading error rate table from ", errorRateFile)
        errorRateTable <- readRDS(errorRateFile)
    }

    if (is.null(errorRateTable))
    {
        errorRateTable <- mutationTable %>%
            filter(SLX_BARCODE %in% patientSamples$SLX_BARCODE) %>%
            mutate(HAS_SIGNAL = ifelse(ALT_F + ALT_R > 0, SLX_BARCODE, NA)) %>%
            group_by(UNIQUE_POS, TRINUCLEOTIDE) %>%
            summarize(MUT_SUM = sum(ALT_F) + sum(ALT_R),
                      DP_SUM = sum(DP),
                      BACKGROUND_AF = MUT_SUM / DP_SUM,
                      N_SAMPLES = n_distinct(SLX_BARCODE),
                      N_SAMPLES_WITH_SIGNAL = n_distinct(HAS_SIGNAL, na.rm = TRUE),
                      .groups = 'drop')

        if (!is.null(errorRateFile))
        {
            saveRDS(errorRateTable, errorRateFile)
        }
    }

    errorRateTable <- errorRateTable %>%
        mutate(PROPORTION = N_SAMPLES_WITH_SIGNAL / N_SAMPLES,
               LOCUS_NOISE.PASS = PROPORTION < proportion_of_controls & BACKGROUND_AF < max_background_mean_AF,
               HAS_AF = BACKGROUND_AF > 0,
               LOCUS_NOISE.FAIL = !LOCUS_NOISE.PASS)

    counts <- errorRateTable %>%
        summarize(NON_ZERO_LOCI_PERCENT = sum(HAS_AF) / n() * 100,
                  LOCUS_NOISE_FAIL_PERCENT = sum(LOCUS_NOISE.FAIL) / n() * 100)

    ## TODO Copy and update plot code.
    ## That said, seems that's for a special case that doesn't apply for most work.

    proportionPositions <- errorRateTable %>%
        filter(PROPORTION < proportion_of_controls)

    backgroundAFPositions <- errorRateTable %>%
        filter(BACKGROUND_AF < max_background_mean_AF)

    extendedMutationTable <- mutationTable %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% proportionPositions$UNIQUE_POS & UNIQUE_POS %in% backgroundAFPositions$UNIQUE_POS)

    extendedMutationTable
}

##
# This function is built from the code after loading in INVAR3.R
#
# take only off_target bases, exclude INDELs
take.offtarget <- function(mutationTable, patientSamples, use_cosmic, is.blood_spot = FALSE)
{
    # Common function. Takes in a mutation table.
    # Adds MUT_SUM, DP_SUM, BACKGROUND_AF columns on rows grouped by
    # SLX, BARCODE, REF, ALT, TRINUCLEOTIDE
    groupAndSummarize <- function(mt)
    {
        mt %>%
            group_by(SLX, BARCODE, REF, ALT, TRINUCLEOTIDE) %>%
            summarize(MUT_SUM = sum(ALT_F) + sum(ALT_R),
                      DP_SUM = sum(DP),
                      BACKGROUND_AF= (sum(ALT_F) + sum(ALT_R)) / sum(DP),
                      .groups = 'drop')

    }

    mutationTable.off_target <- mutationTable %>%
        filter(!ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)

    if (!use_cosmic)
    {
        mutationTable.off_target <- mutationTable.off_target %>%
            filter(!COSMIC)
    }

    # error rate per locus filters
    backgroundErrorTable <- mutationTable.off_target %>%
        groupAndSummarize()

    locusErrorRateFile <- str_c('locus_error_rates.off_target.', ifelse(use_cosmic, 'cosmic', 'no_cosmic'), '.rds')

    mutations.off_target <- annotate_with_locus_error_rate(mutationTable.off_target,
                                                           patientSamples,
                                                           is.blood_spot = is.blood_spot,
                                                           on_target = FALSE,
                                                           errorRateFile = locusErrorRateFile)

    # annotate with both strand information
    mutations.off_target <- mutations.off_target %>%
        mutate(BOTH_STRANDS = (ALT_R > 0 & ALT_F > 0) | AF == 0)

    ## Calculate error rate with differnt settings

    # locus noise pass

    background_error.locus_noise <- mutations.off_target %>%
        filter(LOCUS_NOISE.PASS) %>%
        groupAndSummarize()

    # both strands

    background_error.both_strands <- mutations.off_target %>%
        filter(BOTH_STRANDS) %>%
        groupAndSummarize()

    # locus noise AND both strands

    background_error.locus_noise.both_strands <- mutations.off_target %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS) %>%
        groupAndSummarize()

    ## Assemble all error rates.

    error_rates <- list(pre_filter = backgroundErrorTable,
                        locus_noise_filter_only = background_error.locus_noise,
                        both_strands_only = background_error.both_strands,
                        locus_noise.both_strands = background_error.locus_noise.both_strands)

    error_rates
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationTable <- load.mutations.table(scriptArgs$MUTATIONS_FILE, scriptArgs$MUTATIONS_BED_FILE, scriptArgs$TAPAS_SETTING)

    patientSamples <- load.patient.samples(scriptArgs$LAYOUT_FILE)

    on_target <- filter.for.ontarget(mutationTable)
    saveRDS(on_target, str_c(scriptArgs$TAPAS_SETTING, '.on_target.rds'))
    write_tsv(on_target, str_c(scriptArgs$TAPAS_SETTING, '.on_target.tsv'))

    off_target.cosmic <- take.offtarget(mutationTable, patientSamples, TRUE, scriptArgs$BLOODSPOT)
    saveRDS(off_target.cosmic, str_c(scriptArgs$TAPAS_SETTING, '.off_target.cosmic.error_rates.rds'))

    off_target.no_cosmic <- take.offtarget(mutationTable, patientSamples, FALSE, scriptArgs$BLOODSPOT)
    saveRDS(off_target.no_cosmic, str_c(scriptArgs$TAPAS_SETTING, '.off_target.no_cosmic.error_rates.rds'))
}

# Launch it.

invisible(main(parseOptions()))
