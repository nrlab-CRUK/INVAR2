suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))

source(str_c(Sys.getenv('INVAR_HOME'), '/R/1_parse/common.R'))


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
        make_option(c("--pool"), type="character", metavar="string",
                    dest="POOL", help="The pool (SLX) identifier for the BAM/sizes file",
                    default=defaultMarker),
        make_option(c("--barcode"), type="character", metavar="string",
                    dest="BARCODE", help="The barcode for the BAM/sizes file",
                    default=defaultMarker),
        make_option(c("--outlier-suppression"), type="double", metavar="number",
                    dest="OUTLIER_SUPPRESSION", help="The outlier suppression setting",
                    default=0.05),
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
    list(
        MUTATIONS_TABLE_FILE = 'on_target/mutation_table.on_target.rds',
        FRAGMENT_SIZES_FILE = 'EMMA/output_size/SLX-19721_SXTLI001.inserts_for_annotation.tsv',
        POOL = 'SLX-19721',
        BARCODE = 'SXTLI001',
        OUTLIER_SUPPRESSION = 0.05,
        SAMPLING_SEED = 1024L
    )
}


# Load (locally to Rich) the mutations table from Emma's data into a tibble
# that can be used as if from the pipeline.
readOriginalEmmaMutationsTable <- function(scriptArgs)
{
    combined_polished.path <- '/home/data/INVAR/EMMA/output_R/PARADIGM.f0.9_s2.BQ_20.MQ_40.combined.rds'

    mutationsTable <- as_tibble(readRDS(combined_polished.path)) %>%
        rename(POOL = SLX,
               MUTATION_CLASS = mut_class,
               UNIQUE_POS = uniq_pos) %>%
        filter(POOL == scriptArgs$POOL & BARCODE == scriptArgs$BARCODE)

    mutationsTable
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

        positionSizes <- positionSizes.mutant %>% add_row(positionSizes.wildType)
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

        positionSizes <- positionSizes.mutant %>% add_row(positionSizes.wildType)
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
        add_row(fragmentSizesTable.discrepant.zero.fixed) %>%
        add_row(fragmentSizesTable.discrepant.nonzero.fixed) %>%
        select(UNIQUE_POS, MUTATION_CLASS, SIZE, MUTANT)

    mutationsWithSize <- mutationsTable %>%
        inner_join(fragmentSizesTable.fixed, by = c('UNIQUE_POS', 'MUTATION_CLASS'))

    mutationsWithSize
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        filter(POOL == scriptArgs$POOL & BARCODE == scriptArgs$BARCODE) %>%
        addMutationTableDerivedColumns()

    # mutationsTable <- readOriginalEmmaMutationsTable(scriptArgs)

    fragmentSizesTable <-
        read_tsv(scriptArgs$FRAGMENT_SIZES_FILE, col_types = 'ciccci') %>%
        mutate(MUTANT = ALT == SNV_BASE,
               MUTATION_CLASS = str_c(REF, ALT, sep = '/'),
               UNIQUE_POS = str_c(CHROM, POS, sep = ':'))

    mutationsTable.withSizes <-
        equaliseSizeCounts(mutationsTable, fragmentSizesTable,
                           .samplingSeed = scriptArgs$SAMPLING_SEED)

    mutationsTable.withSizes %>%
        removeMutationTableDerivedColumns() %>%
        arrange(POOL, BARCODE, CHROM, POS, REF, ALT, TRINUCLEOTIDE) %>%
        saveRDSandTSV(str_c('mutation_table.with_sizes.', scriptArgs$POOL, '_', scriptArgs$BARCODE, ".rds"))
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789') {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
