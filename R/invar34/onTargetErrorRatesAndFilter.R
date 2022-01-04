suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggplot2))
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
        make_option(c("--study"), type="character", metavar="string",
                    dest="STUDY", help="The study name",
                    default=defaultMarker),
        make_option(c("--tapas"), type="character", metavar="string",
                    dest="TAPAS_SETTING", help="The TAPAS setting",
                    default=defaultMarker),
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by createMutationsTable.R",
                    default=defaultMarker),
        make_option(c("--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("--cosmic-threshold"), type="integer", metavar="int",
                    dest="COSMIC_THRESHOLD", help="Loci with >0 entries in COSMIC are considered as COSMIC mutations",
                    default=0),
        make_option(c("--bloodspot"), action="store_true", default=FALSE,
                    dest="BLOODSPOT", help="Indicate this is blood spot data."),
        make_option(c("--control-proportion"), type="double", metavar="num",
                    dest="CONTROL_PROPORTION", help="Blacklist loci that have signal in >30% of the nonptspec samples",
                    default=0.1),
        make_option(c("--max-background-af"), type="double", metavar="num",
                    dest="MAX_BACKGROUND_AF", help="Filter loci with a background AF in controls greater than this value",
                    default=0.01),
        make_option(c("--af-threshold"), type="double", metavar="num",
                    dest="AF_THRESHOLD", help="Maximum AF value for acceptable samples",
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
    list(
        MUTATIONS_TABLE_FILE = 'on_target/mutation_table.on_target.all.rds',
        LAYOUT_FILE = 'source_files/combined.SLX_table_with_controls_031220.csv',
        STUDY = 'PARADIGM',
        TAPAS_SETTING = 'f0.9_s2.BQ_20.MQ_40',
        CONTROL_PROPORTION = 0.1,
        MAX_BACKGROUND_AF = 0.01,
        AF_THRESHOLD = 0.01,
        BLOODSPOT = FALSE,
        COSMIC_THRESHOLD = 0
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
        mutate(MUTATION_CLASS = str_c(REF, ALT, sep='/'),
               UNIQUE_POS = str_c(CHROM, POS, sep=':'),
               UNIQUE_ALT = str_c(UNIQUE_POS, MUTATION_CLASS, sep='_'),
               UNIQUE_PATIENT_POS = str_c(PATIENT, UNIQUE_POS, sep='_')) %>%
        select(PATIENT, REF, ALT, TUMOUR_AF, MUTATION_CLASS, UNIQUE_POS, UNIQUE_ALT, UNIQUE_PATIENT_POS)
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
               UNIQUE_ALT = str_c(UNIQUE_POS, str_c(REF, ALT, sep='/'), sep='_'),
               POOL_BARCODE = str_c(POOL, BARCODE, sep='_'),
               MUT_SUM = ALT_F + ALT_R)
}


##
# From TAPAS_functions.R, originally "annotate_with_locus_error_rate"
#
# Originally just annotate, this method creates the lociErrorRateTable and returns.
#
# Locus error rate = overall background error rate per locus, aggregated across control samples
#
createLociErrorRateTable <- function(mutationTable,
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
        filter(!PATIENT_SPECIFIC & CASE_OR_CONTROL == 'case') %>%
        mutate(HAS_SIGNAL = ifelse(ALT_F + ALT_R > 0, POOL_BARCODE, NA)) %>%
        group_by(UNIQUE_POS, MUTATION_CLASS, TRINUCLEOTIDE, PATIENT_MUTATION_BELONGS_TO, COSMIC) %>%
        summarize(MUT_SUM = sum(ALT_F + ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF = MUT_SUM / DP_SUM,
                  N_SAMPLES = n_distinct(POOL_BARCODE),
                  N_SAMPLES_WITH_SIGNAL = n_distinct(HAS_SIGNAL, na.rm = TRUE),
                  .groups = 'drop') %>%
        mutate(LOCUS_NOISE.PASS = (N_SAMPLES_WITH_SIGNAL / N_SAMPLES) < proportion_of_controls &
                                  BACKGROUND_AF < max_background_mean_AF)

    errorRateTable
}


createLociErrorRatePlot <- function(errorRateTable, study, tapasSetting)
{
    nonZeroLoci <- errorRateTable %>%
        summarise(F = sum(BACKGROUND_AF > 0) / n())

    nonZeroLociPercentage <- signif(nonZeroLoci$F[1] * 100, digits = 3)

    locusNoiseFail <- errorRateTable %>%
        summarise(F = sum(!LOCUS_NOISE.PASS) / n())

    locusNoiseFailPercentage <- signif(locusNoiseFail$F[1] * 100, digits = 3)

    plot <- errorRateTable %>%
        ggplot(aes(x = BACKGROUND_AF, fill = LOCUS_NOISE.PASS)) +
        geom_histogram(bins = 100, position = "dodge") +
        scale_colour_discrete(name = "Locus noise pass") +
        scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1)) +
        theme_bw() +
        labs(x = "Background plasma AF across all samples",
             y = "Frequency",
             title = paste(study, tapasSetting, "on-target, non-patient specific data"),
             subtitle = str_c("Blacklisting of non-zero AF loci (representing ", nonZeroLociPercentage,
                              "% of data)\nBlacklisted loci (LOCUS_NOISE.FAIL) = ", locusNoiseFailPercentage,
                              "%\nsplit by COSMIC mutation status")) +
        facet_wrap(~COSMIC, scales = "free_y")

    plot
}


# identify samples with ctDNA >1% so that they are flagged
# later in the pipeline they are not used as patient-controls because of a
# potential risk of contamination between pt-spec and non-pt-spec loci
getContaminatedSamples <- function(mutationTable, afThreshold)
{
    AF <- mutationTable %>%
        mutate(BOTH_STRANDS = ALT_R > 0 & ALT_F > 0 | AF == 0) %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS) %>%
        group_by(SAMPLE_NAME, PATIENT_MUTATION_BELONGS_TO, PATIENT_SPECIFIC) %>%
        summarise(MUT_SUM = sum(ALT_F + ALT_R),
                  DEPTH = sum(DP),
                  AF = MUT_SUM / DEPTH,
                  .groups = 'drop') %>%
        rename(PATIENT = PATIENT_MUTATION_BELONGS_TO)

    ## any correlation between sample and it's nonptspec and ptspec AFs?

    patientSpecific = AF %>%
        filter(PATIENT_SPECIFIC)

    nonPatientSpecific = AF %>%
        filter(!PATIENT_SPECIFIC)

    test <- inner_join(patientSpecific, nonPatientSpecific, by = c('SAMPLE_NAME'))

    lowSamples <- test %>%
        filter(AF.x < afThreshold)

    # print(cor.test(lowSamples$AF.x, lowSamples$AF.y))

    doNotUse <- AF %>%
        filter(PATIENT_SPECIFIC & AF > afThreshold)

    doNotUse %>%
        distinct(SAMPLE_NAME)
}



##
# Saving functions
#

# Remove columns from the mutation table that can be derived from
# other columns, typically before saving.
removeDerivedColumns <- function(mutationTable)
{
    mutationTable %>%
        select(-any_of(c('MUT_SUM', 'POOL_BARCODE')), -contains('UNIQUE'))
}

# Writes the TSV file but, before saving, converts any logical columns to
# simply the characters 'T' or 'F'. Saves having full "TRUE" and "FALSE" values,
# which are excessive as reading the table back correctly interprets 'T' and 'F'.
exportTSV <- function(t, file)
{
    toChar <- function(x)
    {
        ifelse(x, 'T', 'F')
    }

    t %>%
        mutate_if(is.logical, toChar) %>%
        write_tsv(file)

    t
}

# Save the given table as an RDS file and a TSV
saveRDSandTSV <- function(t, file)
{
    saveRDS(t, file)

    tsv <- str_replace(file, "\\.rds$", ".tsv")

    if (tsv != file)
    {
        exportTSV(t, tsv)
    }

    t
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    layoutTable <- loadLayoutFile(scriptArgs$LAYOUT_FILE)

    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addDerivedColumns()

    lociErrorRateTable <- mutationTable %>%
        createLociErrorRateTable(layoutTable,
                                 proportion_of_controls = scriptArgs$CONTROL_PROPORTION,
                                 max_background_mean_AF = scriptArgs$MAX_BACKGROUND_AF,
                                 is.blood_spot = scriptArgs$BLOODSPOT)

    saveRDSandTSV(lociErrorRateTable, 'locus_error_rates.on_target.rds')

    lociErrorRatePlot <- lociErrorRateTable %>%
        createLociErrorRatePlot(study = scriptArgs$STUDY, tapasSetting = scriptArgs$TAPAS_SETTING)

    suppressWarnings(ggsave(lociErrorRatePlot, filename = 'locus_error_rates.on_target.pdf', width = 11, height = 7))

    lociErrorRateNoisePass <- lociErrorRateTable %>%
        filter(LOCUS_NOISE.PASS)

    mutationTable.filtered <- mutationTable %>%
        filter(PATIENT_SPECIFIC | COSMIC_MUTATIONS <= scriptArgs$COSMIC_THRESHOLD) %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% lociErrorRateNoisePass$UNIQUE_POS)

    contaminatedSamples <- getContaminatedSamples(mutationTable.filtered, scriptArgs$AF_THRESHOLD)

    mutationTable.filtered <- mutationTable.filtered %>%
        mutate(CONTAMINATION_RISK.PASS = !SAMPLE_NAME %in% contaminatedSamples$SAMPLE_NAME)

    mutationTable.filtered %>%
        removeDerivedColumns() %>%
        arrange(CHROM, POS, REF, ALT, TRINUCLEOTIDE) %>%
        saveRDSandTSV('mutation_table.on_target.rds')
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789') {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
