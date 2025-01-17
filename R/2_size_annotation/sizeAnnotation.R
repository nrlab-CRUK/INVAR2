suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))


##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"

    options_list <- list(
        make_option(c("--mutations"), type="character", metavar="file",
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) created by part one of the pipeline",
                    default=defaultMarker),
        make_option(c("--fragment-sizes"), type="character", metavar="file",
                    dest="FRAGMENT_SIZES_FILE", help="The fragment sizes file",
                    default=defaultMarker),
        make_option(c("--sample"), type="character", metavar="string",
                    dest="SAMPLE_ID", help="The sample identifier for the BAM/sizes file",
                    default=defaultMarker),
        make_option(c("--outlier-suppression"), type="double", metavar="number",
                    dest="OUTLIER_SUPPRESSION", help="The outlier suppression setting",
                    default=0.05),
        make_option(c("--threads"), type="integer", metavar="integer",
                    dest="THREADS", help="The number of cores to use to process the input files.",
                    default=1L))

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
# Calculation functions.
#

##
# From TAPAS_functions.R "combine_classes_in_rds" but without TRINUCLEOTIDE.
#
# Complement the MUTATION_CLASS column for reference alleles
# 'A' and 'G'. 'T' and 'C' remain unchanged.
#
convertComplementaryMutations <- function(fragmentSizesTable)
{
    complement <- function(sequence)
    {
        chartr('ATCG', 'TAGC', sequence)
    }

    complementary <- function(base)
    {
        length(base) > 0 && (base == 'A' | base == 'G')
    }

    fragmentSizesTable %>%
        rowwise() %>%
        mutate(MUTATION_CLASS = ifelse(complementary(REF), complement(MUTATION_CLASS), MUTATION_CLASS)) %>%
        ungroup()
}

# Resample from a tibble to the given number of rows.
# n == 0 gives an empty table.
# n > rows in 'theTibble' is upsampling: sample with replacement to give more rows.
# n < rows in 'theTibble' is downsampling: sample without replacement to give fewer rows.
# n > 0 when 'theTibble' is empty throws an error.
#
resample <- function(theTibble, n)
{
    assert_that(is.number(n), msg = "n must be a number")
    assert_that(n >= 0, msg = "n must be >= 0")

    tn <- nrow(theTibble)

    sampled <- theTibble

    if (n != tn)
    {
        assert_that(tn > 0, msg = str_c("Requested ", n, " rows from a tibble with no rows."))

        if (n == 0)
        {
            # A way to filter down to an empty tibble keeping the columns.
            sampled <- theTibble %>%
                filter(FALSE)
        }
        else
        {
            # Replace if we want more rows than we have to rows to sample from.
            sampled <- theTibble %>%
                slice_sample(n = n, replace = n > tn)
        }
    }

    sn = nrow(sampled)

    assert_that(sn == n,
                msg = str_c("Have wrong number of rows after sampling: ", sn, " when we wanted ", n))

    sampled
}

# Correct the number of fragment sizes to match the number given
# by mpileup.
#
sampleFragments <- function(pysamFragments, mutationsTableDepths)
{
    grouping <- pysamFragments %>%
        distinct(UNIQUE_POS, MUTANT)

    assert_that(nrow(grouping) == 1,
                msg = "Expect each pysamFragment group to be grouped by UNIQUE_POS and MUTANT")

    uniquePos <- grouping$UNIQUE_POS
    mutant <- grouping$MUTANT
    type <- ifelse(mutant, "mutants", "wild types")

    mpileup <- mutationsTableDepths %>%
        filter(UNIQUE_POS == uniquePos) %>%
        select(all_of(ifelse(mutant, "MUTANT_DEPTH", "WILDTYPE_DEPTH"))) %>%
        rename(DEPTH = 1)

    assert_that(nrow(mpileup) == 1,
                msg = str_c("Expect exactly one row of mpileup depths per unique position: ", uniquePos, ' ', type))

    nfrags <- nrow(pysamFragments)

    if (mpileup$DEPTH > 0 && nfrags == 0)
    {
        stop("Have ", mpileup$DEPTH, " ", type, " from mpileup at ", uniquePos,
             " but have zero rows from pysam.")
    }
    if (mpileup$DEPTH < nfrags)
    {
        message("Have ", mpileup$DEPTH, " ", type, " from mpileup at ", uniquePos,
                " but have ", nfrags, " rows from pysam.")
    }

    pysamFragments.sampled <- pysamFragments %>%
        resample(mpileup$DEPTH)

    nsampled = nrow(pysamFragments.sampled)

    assert_that(nsampled == mpileup$DEPTH,
                msg = str_c("After sampling we don't have the correct number of samples. Expect ",
                            mpileup$DEPTH, " have ", nsampled))

    pysamFragments.sampled
}

# Create a new mutations table that is much bigger where there is a read
# for each matching set of fragments from pysam pileup. There will be quite
# a bit of duplication in the table but each row will have its own MUTANT
# and SIZE value.
#
equaliseSizeCounts <- function(mutationsTable, fragmentSizesTable, threads = 1L)
{
    # Fix depth.

    mutationsTable <- mutationsTable %>%
        mutate(DP = ALT_F + ALT_R + REF_F + REF_R)

    # Summarise with mutant and wild type expected depths

    mutationsTableDepths <- mutationsTable %>%
        mutate(MUTANT_DEPTH = ALT_F + ALT_R,
               WILDTYPE_DEPTH = DP - MUTANT_DEPTH) %>%
        select(UNIQUE_POS, MUTANT_DEPTH, WILDTYPE_DEPTH)

    # Filter for fragments that are in the mutations table only.

    fragmentSizesTable <- fragmentSizesTable %>%
        filter(UNIQUE_POS %in% mutationsTableDepths$UNIQUE_POS)

    # Sample the number of fragments to be the same as the mpileup depth
    # for each position and mutant status.

    fragmentSizesTable.sampled <- fragmentSizesTable %>%
        group_by(UNIQUE_POS, MUTANT) %>%
        group_split(.keep = TRUE) %>%
        mclapply(sampleFragments, mutationsTableDepths, mc.cores = threads) %>%
        bind_rows() %>%
        select(UNIQUE_POS, SIZE, MUTANT)

    # Join with the mutations table to give a table with many more rows
    # (multiply the reads) but each reads will have the SIZE number
    # and MUTANT flag.

    mutationsWithSize <- mutationsTable %>%
        inner_join(fragmentSizesTable.sampled, by = c('UNIQUE_POS'))

    mutationsWithSize
}


##
# Addition to the main function to save the mutation table rows
# for a specific patient (mutation belongs to).
#

saveForPatient <- function(patient, mutationsTable, sampleId)
{
    filename <- str_c('mutation_table.with_sizes', makeSafeForFileName(sampleId), makeSafeForFileName(patient), "rds", sep = '.')

    mutationsTable %>%
        filter(PATIENT_MUTATION_BELONGS_TO == patient) %>%
        removeMutationTableDerivedColumns() %>%
        arrangeMutationTableForExport() %>%
        saveRDS(filename)

    tibble(SAMPLE_ID = sampleId, PATIENT_MUTATION_BELONGS_TO = patient, FILE_NAME = filename)
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    assert_that(file.exists(scriptArgs$MUTATIONS_TABLE_FILE),
                msg = str_c(scriptArgs$MUTATIONS_TABLE_FILE, " does not exist."))

    assert_that(file.exists(scriptArgs$FRAGMENT_SIZES_FILE),
                msg = str_c(scriptArgs$FRAGMENT_SIZES_FILE, " does not exist."))

    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        filter(SAMPLE_ID == scriptArgs$SAMPLE_ID) %>%
        addMutationTableDerivedColumns()

    if (nrow(mutationsTable) == 0)
    {
        warning("There are no rows for sample ", scriptArgs$SAMPLE_ID, " in ",
                basename(scriptArgs$MUTATIONS_TABLE_FILE))
        fileInfoList <- list()
    }
    else
    {
        fragmentSizesTable <-
            read_tsv(scriptArgs$FRAGMENT_SIZES_FILE, col_types = 'ciccci', progress = FALSE) %>%
            mutate(MUTANT = ALT == SNV_BASE,
                   MUTATION_CLASS = str_c(REF, ALT, sep = '/'),
                   UNIQUE_POS = str_c(CHROM, POS, sep = ':'))

        mutationsTable.withSizes <-
            equaliseSizeCounts(mutationsTable, fragmentSizesTable, threads = scriptArgs$THREADS)

        fileInfoList <-
            mclapply(unique(mutationsTable.withSizes$PATIENT_MUTATION_BELONGS_TO), saveForPatient,
                     mutationsTable.withSizes,
                     sampleId = scriptArgs$SAMPLE_ID,
                     mc.cores = scriptArgs$THREADS)
    }

    fileInfoTable <- bind_rows(fileInfoList)

    if (nrow(fileInfoTable) == 0)
    {
        # When there are no rows, we need to write out an empty file info
        # table (with headers). We also need to keep Nextflow happy by writing
        # out a dummy patient file that is otherwise ignored.
        fileInfoTable = tibble(SAMPLE_ID = character(), PATIENT_MUTATION_BELONGS_TO = character(), FILE_NAME = character())

        filename <- str_c('mutation_table.with_sizes', makeSafeForFileName(scriptArgs$SAMPLE_ID), "dummy.rds", sep = '.')

        mutationsTable %>%
            filter(FALSE) %>%
            removeMutationTableDerivedColumns() %>%
            saveRDS(filename)
    }

    if (any(duplicated(fileInfoTable$FILE_NAME)))
    {
        print(fileInfoTable)
        stop("Some file names have clashed after munging. This is a problem that cannot be fixed without changing the patient names.")
    }

    write_csv(fileInfoTable, "index.csv")
}

# Launch it.

invisible(main(parseOptions()))
