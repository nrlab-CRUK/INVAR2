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
                    default=1L),
        make_option(c("--sampling-seed"), type="integer", metavar="number",
                    dest="SAMPLING_SEED", help="The seed for random sampling. Only use in testing.",
                    default=NA))

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
    base <- str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/annotateMutationsWithFragmentSize/source/')

    list(
        MUTATIONS_TABLE_FILE = str_c(base, 'mutation_table.on_target.rds'),
        FRAGMENT_SIZES_FILE = str_c(base, 'SLX-19721.SXTLI001.inserts.tsv'),
        SAMPLE_ID = 'SLX-19721:SXTLI001',
        OUTLIER_SUPPRESSION = 0.05,
        SAMPLING_SEED = 1024L,
        THREADS = 2L
    )
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
        base == 'A' | base == 'G'
    }

    fragmentSizesTable %>%
        mutate(MUTATION_CLASS = ifelse(complementary(REF), complement(MUTATION_CLASS), MUTATION_CLASS))
}

##
# From TAPAS_functions.R
#

# Correct the number of fragment sizes to match the number given
# by mpileup.
#
downsampleFragments <- function(fragmentSizesGroup, uniquePos, .samplingSeed = NA)
{
    assert_that(is.na(.samplingSeed) || is.number(.samplingSeed), msg = ".samplingSeed must be a number or NA.")

    mpileupTotals <- fragmentSizesGroup %>%
        summarise(MUTATION_SUM = ALT_F + ALT_R, DP = DP) %>%
        distinct()

    assert_that(nrow(mpileupTotals) == 1, msg = str_c("Have ", nrow(mpileupTotals), " rows for mpileup totals when there must be exactly 1."))

    positionSizes <- fragmentSizesGroup

    if (summarise(positionSizes, sum(MUTANT)) > mpileupTotals$MUTATION_SUM)
    {
        # discrepant mut_sum, going order the sizes by descending number of entries
        # and set unpaired mutant reads not supported by mpileup to not mutant
        positionSizes.mutant <- positionSizes %>%
            filter(MUTANT) %>%
            arrange(desc(SIZE)) %>%
            mutate(MUTANT = row_number() <= mpileupTotals$MUTATION_SUM)

        # message("Too many mutant position sizes for ", uniquePos, " compared to mpileup: ",
        #         summarise(positionSizes.mutant, n()), " sizes compared to ", mpileupTotals$MUTATION_SUM)
        # message(uniquePos, " set bottom ", summarise(positionSizes.mutant, sum(!MUTANT)), " sizes to be wild type.")

        positionSizes.wildType <- positionSizes %>% filter(!MUTANT)

        positionSizes <-
            bind_rows(positionSizes.mutant, positionSizes.wildType)
    }

    nPositionSizes = summarise(positionSizes, N = n())

    # in some cases, the number of reads from the size CSV differs from the number of mpileup reads
    if (nPositionSizes$N != mpileupTotals$DP)
    {
        positionSizes.mutant <- positionSizes %>% filter(MUTANT)
        positionSizes.wildType <- positionSizes %>% filter(!MUTANT)

        # For exactly reproducible sampling in testing.

        if (!is.na(.samplingSeed))
        {
            set.seed(.samplingSeed)
        }

        # When there are more positions from the python size data than mpileup data
        # we need to downsample wildtypes from the mpileup data (i.e. give fewer rows).

        # When there are fewer positions from the python size data than mpileup data
        # we need to sample with some duplication (replace) to get the number up.
        # Or, in original words, will resample (with replacement) aka stretch mpileup data
        # for wt data (inconsequential as wt reads are not used later in pipeline)

        tooFew = nPositionSizes$N < mpileupTotals$DP

        # message("Too ", ifelse(tooFew, 'few', 'many'),
        #         " position sizes compared to mpileup depth for ", uniquePos,
        #         " - ", nPositionSizes$N, " sizes compared to ", mpileupTotals$DP)

        positionSizes.wildType <- positionSizes.wildType %>%
            slice_sample(n = mpileupTotals$DP - mpileupTotals$MUTATION_SUM, replace = tooFew)

        positionSizes <-
            bind_rows(positionSizes.mutant, positionSizes.wildType)
    }

    positionSizes
}

equaliseSizeCounts <- function(mutationsTable, fragmentSizesTable, .samplingSeed = NA)
{
    assert_that(is.na(.samplingSeed) || is.number(.samplingSeed), msg = ".samplingSeed must be a number or NA.")

    # Fix DP count.

    mutationsTable <- mutationsTable %>%
        mutate(DP = ALT_F + ALT_R + REF_F + REF_R)

    # reduce 12 mutation classes to 6

    fragmentSizesTable <- fragmentSizesTable %>%
        convertComplementaryMutations()

    # summarise the number of mutations and all fragments (DP) at a position.

    fragmentSizesTable.summary <- fragmentSizesTable %>%
        group_by(CHROM, POS, MUTATION_CLASS) %>%
        summarise(MUTATION_SUM.SIZE = sum(MUTANT), DP.SIZE = n(), .groups = 'drop')

    test <- mutationsTable %>%
        left_join(fragmentSizesTable.summary, by = c('CHROM', 'POS', 'MUTATION_CLASS')) %>%
        mutate(CORRECT_DEPTH = DP == DP.SIZE & ALT_F + ALT_R == MUTATION_SUM.SIZE)

    correct <- test %>%
        filter(CORRECT_DEPTH)

    ## correct the raw size df - only AF == 0
    discrepantDP.zero <- test %>%
        filter(!CORRECT_DEPTH & AF == 0)

    # Take the DP top sizes at each locus where AF == 0
    # Also mark all of these positions as not being mutants.
    # This is because the mutant reads identified from size info are not subject to BQ and MQ filters
    fragmentSizesTable.discrepant.zero.fixed <- fragmentSizesTable %>%
        inner_join(select(discrepantDP.zero, UNIQUE_POS, DP), by = "UNIQUE_POS") %>%
        group_by(UNIQUE_POS) %>%
        filter(row_number() <= DP) %>%
        ungroup() %>%
        select(-DP) %>%
        mutate(MUTANT = FALSE)

    ## correct loci with AF > 0
    discrepantDP.nonzero <- test %>%
        filter(!CORRECT_DEPTH & AF > 0)

    fragmentSizesTable.discrepant.nonzero.fixed <- fragmentSizesTable %>%
        inner_join(select(discrepantDP.nonzero, UNIQUE_POS, ALT_F, ALT_R, DP), by = "UNIQUE_POS") %>%
        group_by(UNIQUE_POS) %>%
        group_modify(downsampleFragments, .samplingSeed) %>%
        ungroup() %>%
        select(-ALT_F, -ALT_R, -DP)

    fragmentSizesTable.fixed <- fragmentSizesTable %>%
        filter(UNIQUE_POS %in% correct$UNIQUE_POS) %>%
        bind_rows(fragmentSizesTable.discrepant.zero.fixed) %>%
        bind_rows(fragmentSizesTable.discrepant.nonzero.fixed) %>%
        select(UNIQUE_POS, MUTATION_CLASS, SIZE, MUTANT)

    mutationsWithSize <- mutationsTable %>%
        inner_join(fragmentSizesTable.fixed, by = c('UNIQUE_POS', 'MUTATION_CLASS'))

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
    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        filter(SAMPLE_ID == scriptArgs$SAMPLE_ID) %>%
        addMutationTableDerivedColumns()

    # Expect the combination of pool, barcode and patient mutation belongs to to be
    # unique in this file.

    mutationsFileCheck <- mutationsTable %>%
        distinct(SAMPLE_ID) %>%
        mutate(MATCHING = SAMPLE_ID == scriptArgs$SAMPLE_ID)

    assert_that(nrow(mutationsFileCheck) == 1,
                msg = str_c("Do not have a single pool + barcode in ", scriptArgs$MUTATIONS_TABLE_FILE))

    assert_that(all(mutationsFileCheck$MATCHING),
                msg = str_c("Mutations in", scriptArgs$MUTATIONS_TABLE_FILE, "do not belong to ",
                            scriptArgs$SAMPLE_ID, sep = " "))

    fragmentSizesTable <-
        read_tsv(scriptArgs$FRAGMENT_SIZES_FILE, col_types = 'ciccci') %>%
        mutate(MUTANT = ALT == SNV_BASE,
               MUTATION_CLASS = str_c(REF, ALT, sep = '/'),
               UNIQUE_POS = str_c(CHROM, POS, sep = ':'))

    mutationsTable.withSizes <-
        equaliseSizeCounts(mutationsTable, fragmentSizesTable,
                           .samplingSeed = scriptArgs$SAMPLING_SEED)

    fileInfoList <-
        mclapply(unique(mutationsTable.withSizes$PATIENT_MUTATION_BELONGS_TO), saveForPatient,
                 mutationsTable.withSizes,
                 sampleId = scriptArgs$SAMPLE_ID,
                 mc.cores = scriptArgs$THREADS)

    fileInfoTable <- bind_rows(fileInfoList)

    if (any(duplicated(fileInfoTable$FILE_NAME)))
    {
        print(fileInfoTable)
        stop("Some file names have clashed after munging. This is a problem that cannot be fixed without changing the patient names.")
    }

    write_csv(fileInfoTable, "index.csv")
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
