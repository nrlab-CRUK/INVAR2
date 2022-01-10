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
        make_option(c("--suppression-setting"), type="double", metavar="number",
                    dest="SUPPRESSION_SETTING", help="The suppression setting",
                    default=0.05))

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
        SUPPRESSION_SETTING = 0.05
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

    reverseComplement <- function(sequence)
    {
        stringi::stri_reverse(complement(sequence))
    }

    complementary <- function(base)
    {
        base == 'A' | base == 'G'
    }

    forward <- fragmentSizesTable %>%
        filter(!complementary(REF))

    reverse <- fragmentSizesTable %>%
        filter(complementary(REF))

    reverse <- reverse %>%
        mutate(MUTATION_CLASS = complement(MUTATION_CLASS))

    add_row(forward, reverse) %>%
        arrange(CHROM, POS, REF, ALT)
}

##
# From TAPAS_functions.R
#

downsampleFragments <- function(fragmentSizesGroup)
{
    print(fragmentSizesGroup)

    mpileupTotals <- fragmentSizesGroup %>%
        summarise(MUTATION_SUM = sum(ALT_F + ALT_R),
                  DP = sum(DP))


    positionSizes.mutant <- positionSizes %>% filter(MUTANT)
    positionSizes.wildType <- positionSizes %>% filter(!MUTANT)

    if (nrow(positionSizes.mutant) >= mpileupTotals$MUTATION_SUM)
    {
        # discrepant mut_sum, going order the sizes by descending number of entries
        # and set unpaired mutant reads not supported by mpileup to zero mutant status
        positionSizes.mutant <- positionSizes.mutant %>%
            arrange(desc(SIZE)) %>%
            mutate(MUTANT = MUTANT & row_number() > mpileup.mutationSum)
    }

    positionSizes <- positionSizes.mutant %>% add_rows(positionSizes.wildType)

    # in some cases, the number of reads from the size CSV differs from the number of mpileup reads
    if (nrow(positionSizes) != mpileupTotals$DP)
    {
        positionSizes.mutant <- positionSizes %>% filter(MUTANT)
        positionSizes.wildType <- positionSizes %>% filter(!MUTANT)

        if (nrow(positionSizes) > mpileupTotals$DP)
        {
            # python size data has more rows than mpileup data
            # discrepant DP, downsampling wildtypes from the mpileup data

            positionSizes.wildType <- positionSizes.wildType %>%
                slice_sample(n = mpileup.DP - mpileup.mutationSum, replace = FALSE)
        }
        else
        {
            # python size data has FEWER rows than mpileup data
            # discrepant DP, will resample (with replacement) aka stretch mpileup data
            # for wt data (inconsequential as wt reads are not used later in pipeline)")

            positionSizes.wildType <- positionSizes.wildType %>%
                slice_sample(n = mpileup.DP - mpileup.mutationSum, replace = TRUE)
        }

        positionSizes <- positionSizes.mutant %>% add_rows(positionSizes.wildType)
    }

    positionSizes
}

equaliseSizeCounts <- function(mutationsTable, fragmentSizesTable, pool, barcode)
{
    thisPoolBarcode = str_c(pool, barcode, sep = '_')

    # Filter and fix DP count.

    relevantMutations <- mutationsTable %>%
        filter(POOL_BARCODE == thisPoolBarcode) %>%
        mutate(DP = ALT_F + ALT_R + REF_F + REF_R)

    fragmentSizesTable <- fragmentSizesTable %>%
        mutate(MUTANT = ALT == SNV_BASE,
               MUTATION_CLASS = str_c(REF, ALT, sep = '/'),
               UNIQUE_POS = str_c(CHROM, POS, sep = ':'))

    # reduce 12 mutation classes to 6 and summarise the number of
    # mutations and all fragments (DP) at a position.

    fragmentSizesTable.summary <- fragmentSizesTable %>%
        convertComplementaryMutations() %>%
        group_by(CHROM, POS, MUTATION_CLASS) %>%
        summarise(MUTATION_SUM.SIZE = sum(MUTANT), DP.SIZE = n(), .groups = 'drop')

    test <- relevantMutations %>%
        left_join(fragmentSizesTable.summary, by = c('CHROM', 'POS', 'MUTATION_CLASS')) %>%
        mutate(CORRECT_DEPTH = DP == DP.SIZE & ALT_F + ALT_R == MUTATION_SUM.SIZE,
               UNIQUE_POS = str_c(CHROM, POS, sep = ':'))

    ## correct the raw size df - only AF == 0
    discrepant_DP.zero <- test %>%
        filter(!CORRECT_DEPTH & AF == 0)

    # Take the DP top sizes at each locus where AF == 0
    # Also mark all of these positions as not being mutants.
    # This is because the mutant reads identified from size info are not subject to BQ and MQ filters
    fragmentSizesTable.discrepant.zero.fixed <- fragmentSizesTable %>%
        left_join(select(discrepant_DP.zero, UNIQUE_POS, DP.SIZE), by = "UNIQUE_POS") %>%
        group_by(UNIQUE_POS) %>%
        filter(row_number() <= DP.SIZE) %>%
        ungroup() %>%
        select(-DP.SIZE) %>%
        mutate(MUTANT = FALSE)

    ## correct loci with AF > 0
    discrepant_DP.nonzero <- test %>%
        filter(!CORRECT_DEPTH & AF > 0)

    fragmentSizesTable.discrepant.nonzero.fixed <- fragmentSizesTable %>%
        left_join(select(discrepant_DP.nonzero, UNIQUE_POS, MUTATION_SUM.SIZE, DP.SIZE), by = "UNIQUE_POS") %>%
        group_by(UNIQUE_POS) %>%
        pmap_dfr(downsampleFragments) %>%
        ungroup() %>%
        select(-MUTATION_SUM.SIZE, -DP.SIZE)


    print("1st rbind")
    size_ann.discrepant.fixed <- rbind(size_ann.discrepant.zero.fixed, size_ann.discrepant.nonzero.fixed)

    print("2nd rbind")
    size_ann.fixed <- rbind(size_ann.discrepant.fixed, size_ann[size_ann$uniq_pos %in% correct$uniq_pos,])

    print("merging corrected size df and original combined df")
    size_ann.fixed <- size_ann.fixed %>% dplyr::select(-SLX_barcode)
    combined.size_ann <- inner_join(curr, size_ann.fixed[,3:6], by = c("uniq_pos", "mut_class"))

    return(combined.size_ann)
}


##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns()

    # Mutation status can be R(ef), A(lt) or . (unknown).

    fragmentSizesTable <-
        read_tsv(scriptArgs$FRAGMENT_SIZES_FILE, col_types = 'ciccci')

    equaliseSizeCounts(mutationsTable, fragmentSizesTable,
                       pool = scriptArgs$POOL,
                       barcode = scriptArgs$BARCODE)
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789') {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
