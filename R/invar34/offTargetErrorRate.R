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
        make_option(c("--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
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
        LAYOUT_FILE = 'source_files/combined.SLX_table_with_controls_031220.csv',
        CONTROL_PROPORTION = 0.1,
        MAX_BACKGROUND_AF = 0.01,
        BLOODSPOT = FALSE
    )
}


##
# Loading functions.
#

# Read the layout file and extract unique pool id and barcode pairs.

loadLayoutFile <- function(layoutFile)
{
    suppressWarnings(read_csv(file = layoutFile, col_names = TRUE, show_col_types = FALSE)) %>%
        rename_with(str_to_upper) %>%
        mutate(POOL_BARCODE = str_c(SLX_ID, str_replace(BARCODE, '-', '_'), sep = '_'))
}

addDerivedColumns <- function(mutationTable)
{
    mutationTable %>%
        mutate(UNIQUE_POS = str_c(CHROM, POS, sep=':')) %>%
        mutate(POOL_BARCODE = str_c(POOL, BARCODE, sep='_'))
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

    layoutTable.cases <- layoutTable %>%
        filter(CASE_OR_CONTROL == 'case')

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

# Common function. Takes in a mutation table.
# Adds MUT_SUM, DP_SUM, BACKGROUND_AF columns on rows grouped by
# POOL, BARCODE, REF, ALT, TRINUCLEOTIDE
groupAndSummarizeForErrorRate <- function(mutationTable)
{
    mutationTable %>%
        group_by(POOL, BARCODE, REF, ALT, TRINUCLEOTIDE) %>%
        summarize(MUT_SUM = sum(ALT_F) + sum(ALT_R),
                  DP_SUM = sum(DP),
                  BACKGROUND_AF= (sum(ALT_F) + sum(ALT_R)) / sum(DP),
                  .groups = 'drop')

}

# Filter for rows that are off target, optionally also those that have COSMIC values.

filterForOffTarget <- function(mutationTable, withCosmic)
{
    mutationTable.off_target <- mutationTable %>%
        filter(!ON_TARGET & !SNP & nchar(ALT) == 1 & nchar(REF) == 1)

    if (!withCosmic)
    {
        mutationTable.off_target <- mutationTable.off_target %>%
            filter(!COSMIC)
    }

    mutationTable.off_target
}


# Add locus error rate and both strand information
addLocusNoisePass <- function(mutationTable, errorRateTable)
{
    passed.loci <- errorRateTable %>%
        filter(LOCUS_NOISE.PASS) %>%
        select(UNIQUE_POS)

    mutationTable %>%
        mutate(LOCUS_NOISE.PASS = UNIQUE_POS %in% passed.loci$UNIQUE_POS,
               BOTH_STRANDS = (ALT_R > 0 & ALT_F > 0) | AF == 0)
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

##
# Part of main, where this code is run with and without cosmic.
#

doMain <- function(withCosmic, mutationTable, layoutTable, errorRateTable)
{
    cosmicFilePart = ifelse(withCosmic, 'cosmic', 'no_cosmic')

    mutationTable.off_target <- mutationTable %>%
        filterForOffTarget(withCosmic) %>%
        addLocusNoisePass(errorRateTable)

    all <- mutationTable.off_target %>%
        groupAndSummarizeForErrorRate()

    all %>%
        removeDerivedColums() %>%
        write_tsv(str_c('mutation_table.off_target.', cosmicFilePart, '.oneread.tsv'))

    ## Calculate error rate with different settings and save

    # locus noise pass

    locusNoisePass <- mutationTable.off_target %>%
        filter(LOCUS_NOISE.PASS) %>%
        groupAndSummarizeForErrorRate()

    locusNoisePass %>%
        removeDerivedColums() %>%
        write_tsv(str_c('mutation_table.off_target.', cosmicFilePart, '.locusnoise.tsv'))

    # both strands

    bothStrands <- mutationTable.off_target %>%
        filter(BOTH_STRANDS) %>%
        groupAndSummarizeForErrorRate()

    bothStrands %>%
        removeDerivedColums() %>%
        write_tsv(str_c('mutation_table.off_target.', cosmicFilePart, '.bothreads.tsv'))

    # locus noise AND both strands

    bothFilters <- mutationTable.off_target %>%
        filter(LOCUS_NOISE.PASS & BOTH_STRANDS) %>%
        groupAndSummarizeForErrorRate()

    bothFilters %>%
        removeDerivedColums() %>%
        write_tsv(str_c('mutation_table.off_target.', cosmicFilePart, '.locusnoise_bothreads.tsv'))

    # Save these as a combined RDS file in a list.

    allErrorRates = list(
        ONE_READ = all,
        LOCUS_NOISE = locusNoisePass,
        BOTH_READS = bothStrands,
        LOCUS_NOISE.BOTH_READS = bothFilters
    )

    saveRDS(allErrorRates, str_c('mutation_table.error_rates.', cosmicFilePart, '.rds'))
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    layoutTable <-
        loadLayoutFile(scriptArgs$LAYOUT_FILE) %>%
        filter(CASE_OR_CONTROL == "case")

    mutationTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addDerivedColumns()

    errorRateTable <- createErrorRateTable(mutationTable, layoutTable,
                                           proportion_of_controls = scriptArgs$CONTROL_PROPORTION,
                                           max_background_mean_AF = scriptArgs$MAX_BACKGROUND_AF,
                                           is.blood_spot = scriptArgs$BLOODSPOT)

    saveRDSandTSV(errorRateTable, 'locus_error_rates.off_target.rds')

    #mclapply(c(TRUE, FALSE), doMain, mutationTable, layoutTable, errorRateTable)
    doMain(TRUE, mutationTable, layoutTable, errorRateTable)
    doMain(FALSE, mutationTable, layoutTable, errorRateTable)
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789') {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
