library(readr)

if (!require(invar2tapas))
{
    if (!require(devtools))
    {
        install.packages('devtools', repos='https://cran.ma.imperial.ac.uk')
    }

    # Change to the GitHub location.
    devtools::install_local('/home/bowers01/work/invar/tapas')
    library(invar2tapas)
}


TAPAS.setting <- args[1] #'f0.9_s2.BQ_30.MQ_60'
mutations_bed  <- args[2] #MELR: "~/bed/M2/JXP0128_4_MELR_TAPAS.bed" #AVASTM: A1-7_panels.combined.tidy.180109.bed #LUCID:
SLX_layout_path <- args[3]


load.mutations.table <- function(mutationsFile, patientSpecificFile,
                                 cosmic_threshold = 0, max_DP = 1500, min_ref_DP = 5,
                                 individual_MQSB_threshold = 0.01,
                                 n_alt_alleles_threshold = 3, minor_alt_allele_threshold = 2)
{
    print(paste0("reading in table from ", mutationsFile))

    patient_specific <-
        read_tsv(patientSpecificFile,
                 col_names = c('CHROM', 'START', 'END', 'REF', 'ALT'),
                 col_types = 'ciicc') %>%
        mutate(uniq_pos = str_c(CHROM, END, sep=':'))

    data <- read_tsv(mutationsFile, col_types = 'ciccicddddicciidc')

    # Add columns
    # remove soft-masked repeats, identified by lowercase
    data <- data %>%
        filter(!(str_detect(REF, '[acgt]') | str_detect(ALT, '[acgt]'))) %>%
        mutate(AF = (ALT_F + ALT_R) / DP) %>%
        mutate(SLX_barcode = str_c(SLX, BARCODE, sep='_')) %>%
        mutate(uniq_pos = str_c(CHROM, POS, sep=':')) %>%
        mutate(uniq = str_c(uniq_pos, SLX_barcode, sep='_')) %>%
        mutate(COSMIC = COSMIC_MUTATIONS > cosmic_threshold) %>%
        mutate(SNP = `1KG_AF` > 0) %>%
        mutate(on_target = uniq_pos %in% patient_specific$uniq_pos)


    # save on.target pre-filtering
    data_dt.prefilter <- data %>%
        filter(on_target & nchar(ALT) == 1 & nchar(REF) == 1)

    saveRDS(data_dt.prefilter, file = str_c(TAPAS.setting, ".on_target.prefilter.rds"))

    # Apply filters to blacklist loci
    data <- data %>%
        filter(DP < max_DP & !SNP & REF_R + REF_F >= min_ref_DP & (on_target | !COSMIC))

    data <- blacklist.MQSB(data, individual_MQSB_threshold)
    data <- blacklist.multiallelic(data, n_alt_alleles_threshold, minor_alt_allele_threshold)

    saveRDS(data, file = str_c(TAPAS.setting, '.rds'))

    data
}
