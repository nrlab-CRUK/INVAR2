library(dplyr)
library(readr)
library(stringr)


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

# Filter rows that are on target.

filter.for.ontarget <- function(mutationTable)
{
    mutationTable %>% filter(ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)
}

error.free.positions <- function(mutationTable)
{
    erf <- mutationTable %>% summarise(ERROR_FREE_POSITIONS = (1 - sum(BACKGROUND_AF > 0) / n()) * 100)
    as.numeric(erf)
}

print.error.free.positions <- function(mutationTable)
{
    erf <- error.free.positions(mutationTable)
    message("error-free positions: ", erf)
    erf
}

#' Annotate with locus error rate
#' Locus error rate = overall background error rate per locus, aggregated across control samples
annotate_with_locus_error_rate <- function(mutationTable,
                                           slx_layout_path,
                                           proportion_of_controls = 0.1,
                                           max_background_mean_AF = 0.01,
                                           is.blood_spot = FALSE,
                                           on_target = TRUE, plot = TRUE)
{

    if (is.blood_spot)
    {
        message("sWGS/blood spot mutationTable, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
        max_background_mean_AF <- 1
    }

    patientSamples <-
        read_csv(file = slx_layout_path, col_names = TRUE) %>%
        filter(case_or_control == "case") %>%
        mutate(SLX_BARCODE = str_c(SLX_ID, str_replace(barcode, '-', '_'), sep = '_')) %>%
        select(SLX_BARCODE)

    # This just seems to be a way of cutting down the number of mutations
    if (on_target)
    {
        #message("running code on on-target bases, only using nonptspec dataframe")

        # TODO - What is this "data" column?
        # mutationTable <- filter(data == "nonptspec")

        rFileName <- "locus_error_rates.on_target.rds"
    }
    else
    {
        #message("running code on cluster, using all off-target bases")

        #message("Generated list of all patients across all cohorts that can be used for locus noise filter determination. Total: ", nrow(patientSamples), " cases")

        rFileName <- "locus_error_rates.off_target.rds"
    }

    # To speed up testing.
    if (file.exists(rFileName))
    {
        warning("Reading error rate table from ", rFileName)
        errorRateTable <- readRDS(rFileName)
    }
    else
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

        saveRDS(errorRateTable, rFileName)
    }

    errorRateTable <- errorRateTable %>%
        mutate(PROPORTION = N_SAMPLES_WITH_SIGNAL / N_SAMPLES,
               LOCUS_NOISE.PASS = PROPORTION < proportion_of_controls & BACKGROUND_AF < max_background_mean_AF,
               HAS_AF = BACKGROUND_AF > 0,
               LOCUS_NOISE.FAIL = !LOCUS_NOISE.PASS)

#    nonzero_background_AFs <- errorRateTable %>%
#        filter(LOCUS_NOISE.PASS & BACKGROUND_AF > 0)
#        select(BACKGROUND_AF)

    counts <- errorRateTable %>%
        summarize(NON_ZERO_LOCI_PERCENT = sum(HAS_AF) / n() * 100,
                  LOCUS_NOISE_FAIL_PERCENT = sum(LOCUS_NOISE.FAIL) / n() * 100)

    ## TODO Copy and update plot code.

    proportionPositions <- errorRateTable %>%
        filter(PROPORTION < proportion_of_controls)

    backgroundAFPositions <- errorRateTable %>%
        filter(BACKGROUND_AF < max_background_mean_AF)

    extendedMutationTable <- mutationTable %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% proportionPositions$UNIQUE_POS & UNIQUE_POS %in% backgroundAFPositions$UNIQUE_POS)

    extendedMutationTable
}

#take only off_target bases, exclude INDELs
take.offtarget <- function(mutationTable, slx_layout_path, use_cosmic, is.blood_spot = FALSE)
{
    mutationTable.off_target <- mutationTable %>%
        filter(!ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)

    if (!use_cosmic)
    {
        mutationTable.off_target <- mutationTable.off_target %>%
            filter(!COSMIC)
    }

    # error rate per locus filters
    backgroundErrorTable <- mutationTable.off_target %>%
        group_by(SLX, BARCODE, REF, ALT, TRINUCLEOTIDE) %>%
        summarize(BACKGROUND_AF = (sum(ALT_F) + sum(ALT_R)) / sum(DP))

    # print.error.free.positions(background_error.raw)

    mutations.off_target <- annotate_with_locus_error_rate(mutationTable.off_target,
                                                           slx_layout_path,
                                                           is.blood_spot = is.blood_spot,
                                                           on_target = FALSE)

    mutations.off_target
}

#mutationsFile <- args[1]
#TAPAS.setting <- args[2] #'f0.9_s2.BQ_30.MQ_60'
#mutations_bed  <- args[3] #MELR: "~/bed/M2/JXP0128_4_MELR_TAPAS.bed" #AVASTM: A1-7_panels.combined.tidy.180109.bed #LUCID:
#SLX_layout_path <- args[4]

mutationsFile <- 'mutations/SLX-19722_SXTLI005.mutations.tsv'
mutationsFile <- '/mnt/scratchb/bioinformatics/bowers01/invar_emma/EMMA/output_gz/PARADIGM.f0.9_s2.BQ_20.MQ_40.combined.final.ann.tsv'
TAPAS.setting <- 'f0.9_s2.BQ_30.MQ_60'
mutations_bed <- 'bed/PARADIGM_mutation_list_full_cohort_hg19.bed'
slxLayoutFile <- 'bed/combined.SLX_table_with_controls_031220.csv'

mutationTable <- load.mutations.table(mutationsFile, mutations_bed, TAPAS.setting)

## TODO : REMOVE
# mutationTable <- mutationTable %>% filter(ALT != '.')

on_target <- filter.for.ontarget(mutationTable)

off_target <- take.offtarget(mutationTable, slxLayoutFile, TRUE)
