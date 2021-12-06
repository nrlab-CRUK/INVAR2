library(dplyr)
library(readr)
library(stringr)


##
# Taken from TAPAS_functions.R
# Multiallelic blacklisting.
#

# Filter out rows where the MSQB value is greater than the threshold.

blacklist.MQSB <- function(data, individual_MQSB_threshold)
{
    message("number of unique loci (pre MQSB): ", length(unique(data$UNIQUE_POS)))

    message("applying an indvidiual MQSB threshold of ", individual_MQSB_threshold)

    filtered_data <- data %>%
        filter(MQSB > individual_MQSB_threshold)

    message("number of unique loci (post MQSB): ", length(unique(filtered_data$UNIQUE_POS)))

    filtered_data
}

# Create a list of loci that have more than the given number of alleles.

createMultiallelicBlacklist <- function(data, n_alt_alleles_threshold = 3, minor_alt_allele_threshold = 2)
{
   # Filter for rows with positive AF and
    # determine the number of alt alleles per UNIQUE_POS
    multiallelic <- data %>%
        filter(`1KG_AF` > 0) %>%
        mutate(MUT_SUM = ALT_F + ALT_R) %>%
        select(UNIQUE_POS, ALT, MUT_SUM) %>%
        arrange(UNIQUE_POS)

    # Identify loci with >1 alt alleles (e.g. both A>C and A>T at a position)
    duplicate_positions <- multiallelic %>%
        filter(duplicated(UNIQUE_POS))

    multiallelic <- multiallelic %>%
        filter(UNIQUE_POS %in% duplicate_positions)

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

blacklist.multiallelic <- function(data,
                                   n_alt_alleles_threshold = 3, minor_alt_allele_threshold = 2,
                                   blacklistFile = NULL)
{
    message("starting rows before multiallelic blacklist: ", nrow(data))

    multiallelic_blacklist <- createMultiallelicBlacklist(data, n_alt_alleles_threshold, minor_alt_allele_threshold)

    message("length of multiallelic blacklist is ", nrow(multiallelic_blacklist))

    # save the list of multiallelic loci for reference
    if (!is.null(blacklistFile))
    {
        write_tsv(multiallelic_blacklist, blacklistFile, col_names = FALSE, quote = 'none')
    }

    filtered_data <- data %>%
        filter(!UNIQUE_POS %in% multiallelic_blacklist)

    message("rows after multiallelic blacklist filtering: ", nrow(filtered_data))

    filtered_data
}


##
# Loading functions.
#


# Add columns for SLX + barcode pair, chromosome + position pair,
# and an id that is all of these.

add.identifier.columns <- function(data)
{
    data %>%
        mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':')) %>%
        mutate(SLX_BARCODE = str_c(SLX, BARCODE, sep='_')) %>%
        mutate(UNIQUE_ID = str_c(UNIQUE_POS, SLX_BARCODE, sep='_'))
}


# Load the mutations table from file and filter.

load.mutations.table <- function(mutationsFile, patientSpecificFile, tapasSetting,
                                 cosmic_threshold = 0, max_DP = 1500, min_ref_DP = 5,
                                 individual_MQSB_threshold = 0.01,
                                 n_alt_alleles_threshold = 3, minor_alt_allele_threshold = 2)
{
    message("Reading in table from ", mutationsFile)

    patient_specific <-
        read_tsv(patientSpecificFile,
                 col_names = c('CHROM', 'START', 'END', 'REF', 'ALT'),
                 col_types = 'ciicc') %>%
        mutate(UNIQUE_POS = str_c(CHROM, END, sep=':'))

    data <- read_tsv(mutationsFile, col_types = 'ciccicdddddcciidc')

    # Add SLX & positional columns
    data <- add.identifier.columns(data)

    # Remove soft-masked repeats, identified by lowercase
    # Add columns
    data <- data %>%
        filter(!(str_detect(REF, '[acgt]') | str_detect(ALT, '[acgt]'))) %>%
        mutate(AF = (ALT_F + ALT_R) / DP) %>%
        mutate(COSMIC = COSMIC_MUTATIONS > cosmic_threshold) %>%
        mutate(SNP = `1KG_AF` > 0) %>%
        mutate(ON_TARGET = UNIQUE_POS %in% patient_specific$UNIQUE_POS)


    # save on.target pre-filtering
    data_dt.prefilter <- data %>%
        filter(ON_TARGET & nchar(ALT) == 1 & nchar(REF) == 1)

    saveRDS(data_dt.prefilter, file = str_c(tapasSetting, ".on_target.prefilter.rds"))

    # Apply filters to blacklist loci
    data <- data %>%
        filter(DP < max_DP & !SNP & REF_R + REF_F >= min_ref_DP & (ON_TARGET | !COSMIC))

    data <- blacklist.MQSB(data, individual_MQSB_threshold)
    data <- blacklist.multiallelic(data, n_alt_alleles_threshold, minor_alt_allele_threshold)

    saveRDS(data, file = str_c(tapasSetting, '.rds'))

    data
}

# Filter rows that are on target.

filter.for.ontarget <- function(data)
{
    data %>% filter(ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)
}

error.free.positions <- function(data)
{
    erf <- data %>% summarise(ERROR_FREE_POSITIONS = (1 - sum(BACKGROUND_AF > 0) / n()) * 100)
    as.numeric(erf)
}

print.error.free.positions <- function(data)
{
    erf <- error.free.positions(data)
    message("error-free positions: ", erf)
    erf
}

#' Annotate with locus error rate
#' Locus error rate = overall background error rate per locus, aggregated across control samples
annotate_with_locus_error_rate <- function(data,
                                           slx_layout_path,
                                           proportion_of_controls = 0.1,
                                           max_background_mean_AF = 0.01,
                                           is.blood_spot = FALSE,
                                           on_target = TRUE, plot = TRUE)
{

    if (is.blood_spot)
    {
        message("sWGS/blood spot data, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
        max_background_mean_AF <- 1
    }

    patient_SLX_barcode <-
        read_csv(file = slx_layout_path, col_names = TRUE) %>%
        filter(case_or_control == "case") %>%
        mutate(SLX_BARCODE = str_c(SLX_ID, str_replace(barcode, '-', '_'), sep = '_')) %>%
        select(SLX_BARCODE)

    data <- add.identifier.columns(data)

    if (on_target)
    {
        message("running code on on-target bases, only using nonptspec dataframe")

        # TODO - What is this "data" column?
        # data <- filter(data == "nonptspec")

        rFileName <- "locus_error_rates.on_target.rds"
    }
    else
    {
        message("running code on cluster, using all off-target bases")

        message("Generated list of all patients accross all cohorts that can be used for locus noise filter determination. Total: ", length(patient_SLX_barcode), " cases")

        rFileName <- "locus_error_rates.off_target.Rdata"
    }

    error_rate.raw <- data %>%
       #filter(SLX_BARCODE %in% patient_SLX_barcode) %>%
       mutate(HAS_SIGNAL = ifelse(ALT_F + ALT_R > 0, SLX_BARCODE, NA)) %>%
       summarize(MUT_SUM = sum(ALT_F) + sum(ALT_R),
                 DP = sum(DP),
                 N_SAMPLES = n_distinct(SLX_BARCODE),
                 N_SAMPLES_WITH_SIGNAL = n_distinct(HAS_SIGNAL, na.rm = TRUE)) %>%
       mutate(BACKGROUND_AF = MUT_SUM / DP)

    saveRDS(error_rate.raw, rFileName)

    error_rate.raw <- error_rate.raw %>%
        mutate(PROPORTION = N_SAMPLES_WITH_SIGNAL / N_SAMPLES,
               LOCUS_NOISE.PASS = PROPORTION < proportion_of_controls & BACKGROUND_AF < max_background_mean_AF,
               HAS_AF = BACKGROUND_AF > 0,
               LOCUS_NOISE.FAIL = !LOCUS_NOISE.PASS)

#    nonzero_background_AFs <- error_rate.raw %>%
#        filter(LOCUS_NOISE.PASS & BACKGROUND_AF > 0)
#        select(BACKGROUND_AF)

    counts <- error_rate.raw %>%
        summarize(NON_ZERO_LOCI_PERCENT = sum(HAS_AF) / n() * 100,
                  LOCUS_NOISE_FAIL_PERCENT = sum(LOCUS_NOISE.FAIL) / n() * 100)

    ## TODO Copy and update plot code.

    proportion_positions <- error_rate.raw %>%
        filter(PROPORTION < proportion_of_controls) %>%
        select(UNIQUE_POS)

    background_af_positions <- error_rate.raw %>%
        filter(BACKGROUND_AF < max_background_mean_AF) %>%
        select(UNIQUE_POS)

    data %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% proportion_positions & UNIQUE_POS %in% background_af_positions)
}

#take only off_target bases, exclude INDELs
take.offtarget <- function(data, slx_layout_path, use_cosmic, is.blood_spot = FALSE)
{
    data.off_target <- data %>%
        filter(!ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)

    if (use_cosmic)
    {
        message("NOT removing cosmic loci")
    }
    else
    {
        message("Removing cosmic loci")
        data.off_target <- data.off_target %>%
            filter(!COSMIC)
    }

    # error rate per locus filters
    background_error.raw <- data.off_target %>%
        group_by(SLX, BARCODE, REF, ALT, TRINUCLEOTIDE) %>%
        mutate(MUT_SUM = sum(ALT_F)+sum(ALT_R),
               TOTAL_DP = sum(DP),
               BACKGROUND_AF = MUT_SUM / TOTAL_DP) %>%
        select(SLX, BARCODE, REF, ALT, TRINUCLEOTIDE, BACKGROUND_AF)

    # print.error.free.positions(background_error.raw)

    data_dt.off_target <- annotate_with_locus_error_rate(data.off_target,
                                                         slx_layout_path,
                                                         is.blood_spot = is.blood_spot,
                                                         on_target = FALSE)

    data_dt.off_target
}

#mutationsFile <- args[1]
#TAPAS.setting <- args[2] #'f0.9_s2.BQ_30.MQ_60'
#mutations_bed  <- args[3] #MELR: "~/bed/M2/JXP0128_4_MELR_TAPAS.bed" #AVASTM: A1-7_panels.combined.tidy.180109.bed #LUCID:
#SLX_layout_path <- args[4]

mutationsFile <- 'mutations/SLX-19722_SXTLI005.mutations.tsv'
TAPAS.setting <- 'f0.9_s2.BQ_30.MQ_60'
mutations_bed <- 'bed/PARADIGM_mutation_list_full_cohort_hg19.bed'
slxLayoutFile <- 'bed/combined.SLX_table_with_controls_031220.csv'

data <- load.mutations.table(mutationsFile, mutations_bed, TAPAS.setting)

## TODO : REMOVE
data <- data %>% filter(ALT != '.')

on_target <- filter.for.ontarget(data)

off_target <- take.offtarget(data, slxLayoutFile, TRUE)
