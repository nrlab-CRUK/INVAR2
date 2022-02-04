suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
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
                    dest="MUTATIONS_TABLE_FILE", help="The mutations table file (RDS) with outlier suppression flags",
                    default=defaultMarker),
        make_option(c("--layout"), type="character", metavar="file",
                    dest="LAYOUT_FILE", help="The sequencing layout file",
                    default=defaultMarker),
        make_option(c("--error-rates"), type="character", metavar="file",
                    dest="ERROR_RATES_FILE", help="The off target, no cosmic error rates file",
                    default=defaultMarker),
        make_option(c("--invar-scores"), type="character", metavar="file",
                    dest="INVAR_SCORES_FILE", help="The mutations table file (RDS) created by part one of the pipeline",
                    default=defaultMarker))

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
    testhome <- str_c(Sys.getenv('INVAR_HOME'), '/testing/testdata/')
    base <- str_c(testhome, 'analysis/source/')

    list(
        MUTATIONS_TABLE_FILE = str_c(base, 'mutation_table.with_outliers.rds'),
        INVAR_SCORES_FILE = str_c(base, 'invar_scores.rds'),
        ERROR_RATES_FILE = str_c(base, 'background_error_rates.rds'),
        LAYOUT_FILE = str_c(testhome, 'source_files/combined.SLX_table_with_controls_031220.csv'),
    )
}


cohortMutationContextPlot <- function(contextMutationsTable)
{
    study = unique(contextMutationsTable$STUDY)
    assert_that(length(study) == 1, msg = "Have more than one study in the mutations table.")

    plot <- contextMutationsTable %>%
        ggplot(aes(x = as.character(TRINUCLEOTIDE), fill = MUTATION_CLASS)) +
            geom_bar() +
            theme_classic() +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14,face="bold"),
                  axis.text.x = element_text(angle = 90),
                  legend.position = "bottom") +
            labs(x = "Mutation context",
                 y = "Count",
                 title = str_c("Mutation context in ", study, " cohort"),
                 subtitle = "Unique patient specific mutations in polished filtered data")

    plot
}

cohortMutationAFByClassPlot <- function(contextMutationsTable, layoutTable)
{
    contextMutationsTable <- contextMutationsTable %>%
        filter(TUMOUR_AF > 0)

    study = unique(contextMutationsTable$STUDY)
    assert_that(length(study) == 1, msg = "Have more than one study in the mutations table.")

    plot <- contextMutationsTable %>%
        ggplot(aes(x = TUMOUR_AF, fill = MUTATION_CLASS)) +
            geom_histogram(binwidth = 0.02, position = "stack") +
            theme_classic() +
            labs(x = "AF tumour",
                 y = "Mutation count",
                 title = str_c("Mutation count in ", study, " tumour data"))+
            scale_fill_discrete(name = "Mutation")

    plot
}

mutationsPerPatientPlot <- function(patientSummaryTable)
{
    study = unique(patientSummaryTable$STUDY)
    assert_that(length(study) == 1, msg = "Have more than one study in the mutations table.")

    plot <- patientSummaryTable %>%
        ggplot(aes(x = PATIENT, y = MUTATIONS)) +
            geom_bar(stat = "identity") +
            theme_classic()+
            labs(x = "Patient",
                 y = "Tumour mutations",
                 title = str_c("Tumour mutation count in ", study, " cohort"))

    plot
}

mutationClassByCohortPlot <- function(contextMutationsTable)
{
    contextMutationsClassSummary <- contextMutationsTable %>%
        group_by(STUDY) %>%
        mutate(TOTAL_MUTATIONS = n()) %>%
        group_by(STUDY, MUTATION_CLASS) %>%
        summarise(MUTATIONS = n(),
                  TOTAL_MUTATIONS = unique(TOTAL_MUTATIONS),
                  .groups = "drop") %>%
        mutate(FRACTION = MUTATIONS / TOTAL_MUTATIONS)

    contextMutationsClassSummary %>%
        ggplot(aes(x = mut_class, y = fraction)) +
            geom_bar(stat = "identity")+
            theme_classic() +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14,face="bold"),
                  axis.text.x = element_text(angle = 90),
                  legend.position = c(0.75, 0.8))+
            labs(x = "Mutation class",
                 y = "Fraction",
                 title = str_c("Mutation class in ", study, " panel"),
                 subtitle = "Unique patient specific mutations in polished filtered data")
}


##
# From functions.R, get.error_rate_polishing_comparison
#

errorRatePolishingComparison <- function(errorRatesList, layoutTable, exclude_PPC = T, study , setting, plot = c(F,T)){
    # load background rates
    background_error <- calculate.background_error(error_file, SLX_layout, exclude_PPC = T)
    background_error <- add.missing_error_classes(background_error, trinucleotide_depth)
    background_error <- combine_classes(background_error)

    if (plot == T){
        print("plotting data")
        # set levels
        background_error$data <- factor(background_error$data, levels = c("one_read", "both_reads", "locus_noise", "locus_noise.both_reads"))

        print("plotting individual error rates")
        background_error <<- background_error
        print(filter(background_error, case_or_control == "case", data == "locus_noise.both_reads") %>%
                  ggplot(aes(x = reorder(TRINUCLEOTIDE, background_AF),  y = background_AF, colour = mut_class))+
                  geom_point()+
                  theme_classic()+
                  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
                        axis.text=element_text(size=10),
                        axis.title=element_text(size=14))+
                  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4))+
                  labs(x = "Trinucleotide context",
                       y = "Background error rate",
                       title = paste("Background error rates for", study, setting),
                       subtitle = "Using cases only") +
                  guides(fill=FALSE)+
                  annotation_logticks(sides = "l"))

        ggsave(paste0(plot_dir, "p7_", study, "_", setting, ".", Sys.Date(), ".pdf"), width = 8, height = 4)

        print("Looking into polishing strategies")
        # non changing variables
        error_polishing_settings <- unique(background_error$data)
        one_read_error_rates <- filter(background_error, data == "one_read")
        names(one_read_error_rates)[names(one_read_error_rates) == "background_AF"] <- "background_AF.one_read"

        # looping through the settings
        for (i in 1:length(error_polishing_settings)){
            print(error_polishing_settings[i])

            two_read_error_rates <- filter(background_error, data == error_polishing_settings[i])
            names(two_read_error_rates)[names(two_read_error_rates) == "background_AF"] <- "other_setting"

            error_rates_comparison <- inner_join(one_read_error_rates, two_read_error_rates[,c("TRINUCLEOTIDE", "case_or_control", "ALT", "other_setting")], by = c("TRINUCLEOTIDE", "ALT", "case_or_control"))

            if (error_polishing_settings[i] == "locus_noise.both_reads"){
                print("saving locus_noise.both_reads")
                error_rates_comparison.locus_noise.both_reads <- error_rates_comparison
            }

            error_rates_comparison$ratio <- error_rates_comparison$other_setting / error_rates_comparison$background_AF.one_read

            curr_plot <- ggplot(filter(error_rates_comparison, case_or_control == "case"), aes(x = background_AF.one_read, y = other_setting, color = mut_class))+
                geom_point() +
                scale_y_log10(limits = c(0.5e-7, 1e-2))+
                scale_x_log10(limits = c(0.5e-7, 1e-2))+
                geom_abline(slope = 1, linetype = "dashed")+
                labs(x = "Error rate pre-filter",
                     y = paste("Error rate post-filter"),
                     title = error_polishing_settings[i]) +
                theme_classic()+
                guides(fill=FALSE) +
                theme(axis.text=element_text(size=12),
                      axis.title=element_text(size=14))

            print(curr_plot)

            assign(paste0("error_polishing_", error_polishing_settings[i]),curr_plot)
        }

        print("Now saving all the plots combined")

        ggsave(paste0(plot_dir, "p9_", study, "_", setting,".", Sys.Date(), ".pdf"),
               plot = multiplot(error_polishing_locus_noise,
                                error_polishing_both_reads,
                                #error_polishing_polished_only,
                                error_polishing_locus_noise.both_reads,
                                cols = 1),
               width = 6, height = 9)


        print("Making plot to compare filtering steps in a boxplot - p20")
        #background_error$data <- factor(background_error$data, levels = c("one_read", "both_reads", "polished_only", "polished.both_reads"))
        filter(background_error, case_or_control == "case")%>%
            ggplot(aes(x = mut_class, y = background_AF, fill = data)) +
            geom_boxplot() +
            scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))+
            theme_classic() +
            scale_fill_discrete(name = "Background filtering")+
            labs(x= "Mutation class",
                 y= "Background AF",
                 title= paste("Background error filtering steps\n", study, setting))
        ggsave(paste0(plot_dir, "p20_", study, "_", setting,".",Sys.Date(), ".pdf"), width = 6, height = 4)

    } else{
        print("not plotting data, in order to plot, please set variable in function")
    }
    return(background_error)
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
    layoutTable <-
        loadLayoutTable(scriptArgs$LAYOUT_FILE) %>%
        select(STUDY, POOL_BARCODE, TIMEPOINT)

    mutationsTable <-
        readRDS(scriptArgs$MUTATIONS_TABLE_FILE) %>%
        addMutationTableDerivedColumns() %>%
        mutate(UNIQUE_PATIENT_POS = str_c(PATIENT, UNIQUE_POS, sep="_")) %>%
        filter(PATIENT_SPECIFIC & LOCUS_NOISE.PASS & BOTH_STRANDS & OUTLIER.PASS) %>%
        left_join(layoutTable, by = c('STUDY', 'POOL_BARCODE'))

    invarScoresTable <- readRDS(scriptArgs$INVAR_SCORES_FILE)


    contextMutationsTable <- mutationsTable %>%
        group_by(PATIENT, UNIQUE_POS) %>%
        slice_head(n = 1) %>%
        ungroup()

    patientSummaryTable <- contextMutationsTable %>%
        group_by(STUDY, PATIENT) %>%
        summarise(MUTATIONS = n_distinct(STUDY, PATIENT, UNIQUE_POS), .groups = "drop") %>%
        arrange(PATIENT, MUTATIONS)

    errorRatesList <- readRDS(scriptArgs$ERROR_RATES_FILE)

    patientSummaryTable %>%
        write_csv("tumour_mutation_per_patient.csv")

    # plotting 3bp context

    ggsave(plot = cohortMutationContextPlot(contextMutationsTable),
           filename = "p1_cohort_mut_context.pdf",
           width = 6, height = 5)

    # plot the AF by mutation class

    ggsave(plot = cohortMutationAFByClassPlot(contextMutationsTable),
           filename = "p2_cohort_mut_AF_by_class.pdf",
           width = 6, height = 7)

    # plot mutations per patient captured and passing pipeline filters

    ggsave(plot = mutationsPerPatientPlot(patientSummaryTable),
           filename = "p3_cohort_mut_count_tumour.pdf",
           width = 6, height = 7)

    # plot mutation class distribution by cohort

    ggsave(plot = mutationClassByCohortPlot(contextMutationsTable),
           filename = "p4_cohort_mut_class.pdf",
           width = 6, height = 5)

    stop("Stop here")


    #You will have to reload the settings section here
    #{r background error rate}
    # get path to the file of interest:
    file <- list.files(input_dir, pattern = "cosmic_FALSE.error_rates.Rdata", full.names = T)

    # get summary plots
    background_error <- get.error_rate_polishing_comparison(file, SLX_layout, study = study, setting = "fam_2", plot = T)

    # get error rates by class by case vs. control
    load(file)

    error_rates <- error_rates[["pre_filter"]]
    error_rates.parsed <- calculate_background_error_INV042(error_rates = error_rates, SLX_layout = SLX_layout)
    error_rates.parsed$mut_class <- paste0(error_rates.parsed$REF, "/", error_rates.parsed$ALT)

    error_rates.complement1 <- filter(error_rates.parsed, REF == "C" | REF== "T")
    error_rates.complement2 <- filter(error_rates.parsed, REF == "A" | REF== "G")

    error_rates.complement2$mut_class <- chartr("ATGC","TACG", error_rates.complement2$mut_class)

    error_rates.parsed <- rbind(error_rates.complement1, error_rates.complement2)

    ggplot(error_rates.parsed, aes(x = mut_class, y= background_AF, colour = case_or_control))+
      geom_boxplot()+
      scale_y_log10()+
      stat_compare_means(method = "t.test",
                        method.args = list(alternative = "greater"),
                        aes(label = paste0("p = ", ..p.format..)))+
      theme_bw()+
      labs(x = "Mutation class",
           y = "Background error rate",
           title = "Comparison of error rate estimates in case vs. control",
           subtitle = paste(study, "cohort"))+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14),
            panel.grid.major = element_line(colour = alpha("black", 0.1)))+
      stat_compare_means(method = "t.test",
                         method.args = list(alternative = "greater"),
                         aes(label = paste0("p = ", ..p.format..)))
    ggsave(paste0(plot_dir, "p5_background_error_comparison.case_vs_control.pdf"), width = 6, height = 4)

    ## min and max error rate by class
    error_rates.by_class <- plyr::ddply(error_rates.parsed, "mut_class", function(x){
      data.frame(error_rate = weighted.mean(x$background_AF))
    })

    min(error_rates.by_class$error_rate)
    max(error_rates.by_class$error_rate)

    ## min and max error rate by trinucleotide context
    error_rates.by_context <- plyr::ddply(error_rates.parsed, c("mut_class", "TRINUCLEOTIDE"), function(x){
      data.frame(error_rate = weighted.mean(x$background_AF))
    })
    min(error_rates.by_context$error_rate)
    max(error_rates.by_context$error_rate)

    error_rate.summary <- data.frame(setting = c("min_by_class", "max_by_class", "min_by_3bp", "max_by_3bp"),
                                     error_rate = c(min(error_rates.by_class$error_rate), max(error_rates.by_class$error_rate),
                                                    min(error_rates.by_context$error_rate), max(error_rates.by_context$error_rate)))

    write.csv(error_rate.summary, paste0(output_dir, "error_rate.summary.csv"), row.names = F)


    #{r outlier-suppression}
    files <- list.files(input_dir, pattern = "f0.9_s2.BQ_20.MQ_40.combined.os_0.05.rds", full.names = T)

    summary.cohorts <- data.frame()
    output.raw <- data.frame()
    cont_risk_samples <- data.frame()

    for (i in 1:length(files)) {
      print(files[i])
      curr <- readRDS(files[i])
      # samples with high ctDNA levels will not be used for the outlier suppression analysis and won't be plotted in here. These are:
      print(paste0("contamination risk samples that won't be plotted: ", unique(filter(curr, CONTAMINATION_RISK.PASS == FALSE)$SLX_barcode)))
      cont_risk_samples.curr <- as.data.frame(unique(filter(curr, CONTAMINATION_RISK.PASS == FALSE)$SLX_barcode))
      cont_risk_samples <- rbind(cont_risk_samples, cont_risk_samples.curr)
      curr.filtered <- filter(curr,
                              LOCUS_NOISE.PASS == TRUE,
                              BOTH_STRANDS == TRUE,
                              CONTAMINATION_RISK.PASS == TRUE,
                              AF> 0,
                              AF < 0.25,
                              mut_sum < 10)

      # rbind the raw output for the plot of individual data points
      output.raw <- rbind(curr.filtered, output.raw)

      ## get proportion of nonzero loci that were retained by PASS filter
      summary.cohorts.curr <- plyr::ddply(curr.filtered, "data", function(x){
         nrow_pass <- nrow(filter(x, PASS == TRUE))
         data.frame(proportion = nrow_pass/nrow(x),
                    total_rows = nrow(x),
                    study = unique(x$Study))
      })

      summary.cohorts <- rbind(summary.cohorts,
                               summary.cohorts.curr)
    }

    mean(filter(summary.cohorts, data == "ptspec")$proportion)
    mean(filter(summary.cohorts, data == "nonptspec")$proportion)
    summary.cohorts[summary.cohorts$data == "nonptspec",'data'] <- "controls"

    summary.cohorts$data <- factor(summary.cohorts$data, levels = c("ptspec", "controls"))

    ggplot(summary.cohorts, aes(x = data, y = proportion, fill = study, label = round(proportion, digits = 2)))+
      geom_bar(position = "dodge", stat = "identity")+
      geom_text(vjust = 1.2)+
      labs(title = paste0(study, ": Effect of outlier-suppression"),
           x = "Data",
           y = "Proportion of signal\nretained after filter")+
      theme_classic()+
      theme(legend.position = c(0.7, 0.7))+
      scale_y_continuous(limits = c(0,1))
    ggsave(paste0(plot_dir, "p6_os.data_retained.threshold_0.05.pdf"), width = 5, height = 6)

    ggplot(output.raw, aes(x = sample_name, y = AF, colour = PASS))+
      geom_point()+
      theme_classic()+
      scale_y_log10()+
      facet_wrap( Study ~data, scales = "free_x", nrow = 3)+
      theme(axis.text.x=element_text(angle = 90))+
      labs(title = paste("Effect of outlier-suppression", study),
           subtitle = paste("f0.9_s2, threshold = 0.05"),
           x = "Sample",
           y = "Mutant allele fraction")
    ggsave(paste0(plot_dir, "p8_background_polishing", ".threshold_0.05.",study, ".", setting, ".pdf"), width = 15, height = 12)
    #rm(list=ls())


    #You will have to reload the settings section here
    #{r tumour vs plasma AFs}
    ##### need zero AF data as well...
    rds_files <- list.files(input_dir, pattern = "f0.9_s2.BQ_20.MQ_40.combined.os_0.05.rds", full.names = T)

    combined_data <- data.frame()
    for (k in 1:length(rds_files)) {
      curr.file <- readRDS(rds_files[k])

      curr.file <- filter(curr.file,
           PASS == TRUE,
           BOTH_STRANDS == TRUE,
           CONTAMINATION_RISK.PASS == TRUE,
           LOCUS_NOISE.PASS == TRUE)

      combined_data <- rbind(combined_data, curr.file)
    }

    combined_data$observed_in_plasma <- combined_data$AF > 0

    # recalculate average to exclude high level samples
    average <- filter(combined_data, Study == study)%>%
      plyr::ddply(c("Patient", "SLX_barcode"), function(x){
      data.frame(average = weighted.mean(x$AF, x$DP))
    })

    to_remove <- average[average$average > 1e-2,'SLX_barcode']
    output_filtered <- filter(combined_data, !(SLX_barcode %in% to_remove))

    ggplot(data = output_filtered, aes(x = data, y = tumour_AF, fill = observed_in_plasma))+
      geom_boxplot(outlier.colour = NA, notch = T)+
      theme_classic()+
      labs(x = "Data",
           y = "Tumour allele fraction",
           title = paste0(study, ": Tumour AF of detected vs. non detected loci"),
           subtitle = "Split by ptspec and nonptspec")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            panel.grid.major = element_line(colour = alpha("black", 0.1)))+
      stat_compare_means(method = "t.test",
                         method.args = list(alternative = "greater"),
                         aes(label = paste0("p = ", ..p.format..)))
    ggsave(paste0(plot_dir, "p10_tumour_AF_for_observed_non_observed_loci.pdf"), width = 5, height = 4)


    #{r Size profiles & enrichment ratios}
    # Do size comparison for consensus threshold 0.9 of the different family sizes.
    size_files <- list.files(input_dir, full.names = T, pattern = "BQ_20.MQ_40.size_characterisation.rds")

    summary <- size_characterisation(size_files, round = 5, combine_AVASTM = F)
    summary$family_size <- setting
    summary$study_name <- study

    ## fragment size per cohort, with different levels of error-suppression
    filter(summary, data == "ptspec") %>%
    ggplot( aes(x = size.round, y = proportion, fill = mut)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Fragment size in bp",
           y = "Proportion",
           title = paste0(study, ": Fragment size vs. error-suppression")) +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14,
                                      face = "bold"),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12,
                                        face = "bold")) +
      scale_x_continuous(limits = c(0, 500))+
      facet_grid(study_name~., scales = "free_y")+
      geom_vline(xintercept = c(166, 166*2), linetype = "dashed", alpha = 0.5)+
      annotate("text", x = 200, y = 0.15, label = "166bp", alpha = 0.5)+
      annotate("text", x = 370, y = 0.15, label = "332bp", alpha = 0.5)
    ggsave(paste0(plot_dir, "p11_size_comparison.pdf"), width = 6, height = 5)

    ## calculate enrichment level
    size.mutant  <- filter(summary, mut == TRUE, exome == FALSE)[,c("count", "total", "size.round", "data", "case_or_control", "sample_type", "proportion", "family_size", "study_name")]
    size.wt  <- filter(summary, mut == FALSE, exome == FALSE)[,c("count", "total", "size.round", "data", "case_or_control", "sample_type", "proportion", "family_size", "study_name")]
    enrichment_ratio <- left_join(size.wt, size.mutant, by = c("size.round", "data", "case_or_control", "sample_type", "family_size", "study_name"))

    # enrichment ratio correction for comparing different sized dataframes
    totals <- enrichment_ratio[!is.na(enrichment_ratio$proportion.y),] %>%
      plyr::ddply(c("study_name", "data", "case_or_control", "family_size"), function(x){
        data.frame(total = sum(x$count.x)) })

    enrichment_ratio <- left_join(enrichment_ratio, totals, by = c("study_name", "data", "case_or_control", "family_size"))
    enrichment_ratio$proportion.x <- enrichment_ratio$count.x / enrichment_ratio$total

    enrichment_ratio$ratio <- enrichment_ratio$proportion.y/enrichment_ratio$proportion.x
    enrichment_ratio$ratio.log2 <- log2(enrichment_ratio$ratio)

    filter(enrichment_ratio, data == "ptspec", case_or_control == "case", ratio.log2 != "NA")%>%
      ggplot( aes(x = size.round, y = ratio.log2)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Fragment size in bp",
           y = "Enrichment ratio - log2 ratio",
           title = paste0(study, ": Enrichment ratio for ctDNA vs. error-suppression")) +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14,
                                      face = "bold"))+
      scale_x_continuous(limits = c(0, 500))+
      facet_grid(study_name~.)+
      geom_vline(xintercept = c(166, 166*2), linetype = "dashed", alpha = 0.5)+
      annotate("text", x = 210, y = 5, label = "166bp", alpha = 0.5)+
      annotate("text", x = 380, y = 5, label = "332bp", alpha = 0.5)
    ggsave(paste0(plot_dir, "p12_enrichment_ratios.pdf"), width = 6, height = 5)
    #rm(list=ls())


    #You will have to reload the settings section here
    #{r format GLRT and plot ROC}
    INVAR_score.files <<- list.files(input_dir, full.names = T, pattern = "f0.9_s2.BQ_20.MQ_40.INVAR_scores")

    combined_ptspec_GLRT <- data.frame()

    for(x in 1:length(INVAR_score.files)){
      CONTAMINATION_RISK_THRESHOLD = 0.01
      print(paste("loading", INVAR_score.files[x]))
      INVAR_score <- get.INVAR_score(INVAR_score.files[x],
                                     SLX_layout,
                                     adjust = T,
                                     outlier_suppression = 0.05,
                                     filter.pattern = "TRUE_TRUE_TRUE")

      # determine ROC threshold when only using cases, dim(number of bam input, 14)
      ptspec_GLRT <- scale_INVAR_size(INVAR_score, SLX_layout = SLX_layout, downsampled = F)[[1]]
      ptspec_GLRT$data <- INVAR_score.files[x]

      ROC_pt_control <- scale_INVAR_size(INVAR_score, SLX_layout = SLX_layout, downsampled = F)[[3]]

      # names is approx num2str()
      Specificity_pt_control <- names(INVAR_SCORE_threshold.size)
      pt_cutoff <- INVAR_SCORE_threshold.size ## INVAR cut-off value

      # define ROC threshold for using healthy controls if these are present
      INVAR_score_healthy <- filter(INVAR_score, data == "ptspec" | case_or_control == "control_negative")

      if ((nrow(filter(INVAR_score_healthy, case_or_control == "control_negative")) > 0) == TRUE){
        print("this data has healthy controls. Continue with analysis")

        INVAR_score_healthy$case_or_control <- "case"

        ROC_healthy_control <- scale_INVAR_size(INVAR_score_healthy, SLX_layout = SLX_layout, downsampled = F)[[3]]
        Specificity_healthy_control <- names(INVAR_SCORE_threshold.size) # quantile value from scale_INVAR_size()
        healthy_cutoff <- INVAR_SCORE_threshold.size ## INVAR cut-off value

          # plot the combined ROC curve
        ggplot(ROC_pt_control,aes(d = status, m = adjusted_INVAR.y)) +
          geom_roc(data = ROC_healthy_control, cutoffs.at = healthy_cutoff, labels = FALSE, pointsize = 0.75, color = "red")+ # to get rid of numbers on the ROC curve
          geom_roc(cutoffs.at = pt_cutoff, labels = FALSE, pointsize = 0.75)+ # to get rid of numbers on the ROC curve
          theme_classic()+
          labs(title = paste0(study, ", ", setting),
              x = "False positive fraction",
              y = "True positive fraction",
              subtitle = paste0("Specificity patients = ", Specificity_pt_control, "\nSpecificity healthy = ", Specificity_healthy_control))
      ggsave(paste0(plot_dir, "p13_ROC_", study, "_", x, ".pdf"), width = 4, height = 3)

      } else {
        print("This data does not have healthy controls, will compute the ROC curve from case data only")

          # plot the combined ROC curve
        ggplot(ROC_pt_control,aes(d = status, m = adjusted_INVAR.y)) +
          geom_roc(cutoffs.at = pt_cutoff, labels = FALSE, pointsize = 0.75)+ # to get rid of numbers on the ROC curve
          theme_classic()+
          labs(title = paste0(study, ", ", setting),
              x = "False positive fraction",
              y = "True positive fraction",
              subtitle = paste0("Specificity patients = ", Specificity_pt_control))
        ggsave(paste0(plot_dir, "p13_ROC_", study, "_", x, ".pdf"), width = 4, height = 3)

      }

      combined_ptspec_GLRT <- rbind(combined_ptspec_GLRT, ptspec_GLRT)
    }

    #You will have to reload the settings section here
    #{r format GLRT and plot ROC} NO SIZE DATA
    INVAR_score.files <<- list.files(input_dir, full.names = T, pattern = "f0.9_s2.BQ_20.MQ_40.INVAR_scores")

    combined_ptspec_GLRT <- data.frame()

    for(x in 1:length(INVAR_score.files)){
      CONTAMINATION_RISK_THRESHOLD = 0.01
      print(paste("loading", INVAR_score.files[x]))
      INVAR_score <- get.INVAR_score(INVAR_score.files[x],
                                     SLX_layout,
                                     adjust = T,
                                     outlier_suppression = 0.05,
                                     filter.pattern = "TRUE_TRUE_TRUE")

      # determine ROC threshold when only using cases, dim(number of bam input, 14)
      ptspec_GLRT <- scale_INVAR_size(INVAR_score, SLX_layout = SLX_layout, downsampled = F)[[1]]
      ptspec_GLRT$data <- INVAR_score.files[x]

      ROC_pt_control_ns <- scale_INVAR_size(INVAR_score, SLX_layout = SLX_layout, downsampled = F)[[4]]

      # names is approx num2str()
      Specificity_pt_control_ns <- names(INVAR_SCORE_threshold.no_size)
      pt_cutoff_ns <- INVAR_SCORE_threshold.no_size ## INVAR cut-off value

      # define ROC threshold for using healthy controls if these are present
      INVAR_score_healthy <- filter(INVAR_score, data == "ptspec" | case_or_control == "control_negative")

      if ((nrow(filter(INVAR_score_healthy, case_or_control == "control_negative")) > 0) == TRUE){
        print("this data has healthy controls. Continue with analysis")

        INVAR_score_healthy$case_or_control <- "case"

        ROC_healthy_control_ns <- scale_INVAR_size(INVAR_score_healthy, SLX_layout = SLX_layout, downsampled = F)[[4]]
        Specificity_healthy_control_ns <- names(INVAR_SCORE_threshold.no_size) # quantile value from scale_INVAR_size()
        healthy_cutoff_ns <- INVAR_SCORE_threshold.no_size ## INVAR cut-off value

        # plot the combined ROC curve
        ggplot(ROC_pt_control_ns,aes(d = status, m = adjusted_INVAR.x)) +
          geom_roc(data = ROC_healthy_control_ns, cutoffs.at = healthy_cutoff_ns, labels = FALSE, pointsize = 0.75, color = "red")+ # to get rid of numbers on the ROC curve
          geom_roc(cutoffs.at = pt_cutoff_ns, labels = FALSE, pointsize = 0.75)+ # to get rid of numbers on the ROC curve
          theme_classic()+
          labs(title = paste0(study, ", ", setting),
               x = "False positive fraction",
               y = "True positive fraction",
               subtitle = paste0("Specificity patients = ", Specificity_pt_control_ns, "\nSpecificity healthy = ", Specificity_healthy_control_ns))
        ggsave(paste0(plot_dir, "p13a_ROC_", study, "_", x, "no_size.pdf"), width = 4, height = 3)

      } else {
        print("This data does not have healthy controls, will compute the ROC curve from case data only")

        # plot the combined ROC curve
        ggplot(ROC_pt_control_ns,aes(d = status, m = adjusted_INVAR.x)) +
          geom_roc(cutoffs.at = pt_cutoff_ns, labels = FALSE, pointsize = 0.75)+ # to get rid of numbers on the ROC curve
          theme_classic()+
          labs(title = paste0(study, ", ", setting),
               x = "False positive fraction",
               y = "True positive fraction",
               subtitle = paste0("Specificity patients = ", Specificity_pt_control_ns))
        ggsave(paste0(plot_dir, "p13a_ROC_", study, "_", x, "no_size_data.pdf"), width = 4, height = 3)

      }

      combined_ptspec_GLRT <- rbind(combined_ptspec_GLRT, ptspec_GLRT)
    }


    ##### get number of original mutations per patient in the panel, use the newly generated table from upper chunk:
    # will use the mutations that passed the Agilent stage for normalisation
    mutation_list_output <- read.csv(paste0(output_dir, "tumour_mutation_per_pt.csv"), header = T)
    ## get pt id
    combined_ptspec_GLRT$patient <- unlist(lapply(strsplit(combined_ptspec_GLRT$sample_name, split = " "), "[[", 1))

    IF_patient_data <- left_join(combined_ptspec_GLRT, mutation_list_output, by = ("patient" = "Patient"))
    IF_patient_data <- left_join(IF_patient_data, SLX_layout[,c("SLX_barcode", "Study")], by = "SLX_barcode")
    IF_patient_data$uniq_mol <- IF_patient_data$DP/IF_patient_data$mutations
    IF_patient_data$ng_on_seq <- IF_patient_data$uniq_mol/300
    IF_patient_data$low_sensitivity <- IF_patient_data$DP < 20000 & IF_patient_data$detected_using_size == FALSE

    ## stats for paper
    total_DP_per_pt_timepoint <- plyr::ddply(IF_patient_data, c("patient", "Timepoint"), function(x){
      data.frame(total_DP = sum(x$DP))
    })

    median(filter(total_DP_per_pt_timepoint, total_DP > 20000)$total_DP)

    ### effect of filtering
    threshold_settings <- c(0,20000,66666)
    final_output <- data.frame()
    for (k in 1:length(threshold_settings)) {

      curr_detected <- filter(IF_patient_data, DP > threshold_settings[k],
                              detected_using_size == TRUE)%>%
        nrow()

      curr_non_detected <- filter(IF_patient_data, DP > threshold_settings[k],
                              detected_using_size == FALSE)%>%
        nrow()

      curr_cases <- filter(IF_patient_data, DP > threshold_settings[k])%>%
        nrow()

      rm_non_detected <- filter(IF_patient_data, DP < threshold_settings[k],
                            detected_using_size == FALSE)%>%
      nrow()

      output <- data.frame(threshold = threshold_settings[k], curr_cases, curr_detected, curr_non_detected, rm_non_detected)
      final_output <- rbind(final_output, output)
    }

    final_output$all_cases <- nrow(IF_patient_data)
    final_output$all_detected <- filter(IF_patient_data, detected_using_size == TRUE)%>%
        nrow()

    final_output$detection_rate_harsh <- final_output$curr_detected / final_output$curr_cases
    final_output$detection_rate_soft <- final_output$all_detected / (final_output$all_detected + final_output$curr_non_detected)

    write.csv(combined_ptspec_GLRT, paste0(output_dir, "combined_ptspec_GLRT.csv"), row.names = F)
    write.csv(IF_patient_data, paste0(output_dir, "IF_patient_data.csv"), row.names = F)
    write.csv(final_output, paste0(output_dir, "IR_threshold_effects.csv"), row.names = F)


    #{r plotting with GLRT values}
    ############## plot IR (DP) vs IMAF
    IF_patient_data <- read.csv(paste0(output_dir, "IF_patient_data.csv"), header = T)

    IF_patient_data$logable_INVAR <- IF_patient_data$adjusted_IMAF
    IF_patient_data[IF_patient_data$IMAF == 0, "logable_INVAR"] <- 1e-7

    tmp_20k <- data.frame(x = c(5e3,5e3,20000,20000),
                   y = c(0, 2e-04, 5e-05, 0))

    tmp_66k <- data.frame(x = c(20000,20000, 66666, 66666),
                   y = c(0,5e-05, 1.5e-5, 0))

    tmp_1e6 <- data.frame(x = c(1e6,1e6, 5e6, 5e6),
                   y = c(1e-06, 3e-1, 3e-1, 2e-07))

    DP_v <- data.frame(c(5e3, 5e4, 5e5, 5e6))
    colnames(DP_v) <- c("DP")
    IMAF_v <- 1/DP_v
    colnames(IMAF_v) <- c("IMAF")
    sens_line <- cbind(DP_v, IMAF_v)

    ggplot(IF_patient_data, aes(x = DP, y = logable_INVAR))+
      geom_polygon(data = tmp_20k, aes(x = x, y = y), alpha = 0.2, fill = "darkslategray4")+
      geom_polygon(data = tmp_66k, aes(x = x, y = y), alpha = 0.2, fill = "darkslategray3")+
      geom_polygon(data = tmp_1e6, aes(x = x, y = y), alpha = 0.3, fill = "#F59C20")+
      geom_vline(xintercept = 20000, linetype = "dotted")+
      geom_vline(xintercept = 66666, linetype = "dotted")+
      geom_vline(xintercept = 1e6, linetype = "dotted")+
      geom_point(aes(color = detected_using_size))+
      geom_line(data = sens_line, aes(x = DP, y = IMAF),
                  linetype = "longdash")+
      scale_x_log10(limits = c(5e3, 5e6))+
      scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
                    labels = c("ND", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1"))+
      theme_classic()+
      labs(x = "IR",
           y = "IMAF",
           title = "PBCP IR vs IMAF plot")
    ggsave(paste0(plot_dir, "p14_IR_vs_IMAF.pdf"), width = 6, height = 4)

    ############ sample waterfall plot
    combined_ptspec_GLRT <- read.csv(paste0(output_dir, "combined_ptspec_GLRT.csv"), header = T)

    ######### waterfall plot
    combined_ptspec_GLRT$study <- study

    ## add detectable vs non detectable using dPCR
    combined_ptspec_GLRT.annotated <- left_join(combined_ptspec_GLRT, SLX_layout[,c("input_into_library_ng", "quantification_method", "SLX_barcode", "Study")], by = "SLX_barcode")

    combined_ptspec_GLRT.annotated$not_detectable_dPCR <- combined_ptspec_GLRT.annotated$adjusted_IMAF < 3.3/combined_ptspec_GLRT.annotated$input_into_library_ng

    combined_ptspec_GLRT.annotated$ctDNA_plotting <- combined_ptspec_GLRT.annotated$adjusted_IMAF
    combined_ptspec_GLRT.annotated[combined_ptspec_GLRT.annotated$detected_using_size == F, "ctDNA_plotting"] <- 1e-07

    combined_ptspec_GLRT.annotated$ctDNA_plotting.log <- log10(combined_ptspec_GLRT.annotated$ctDNA_plotting)

    combined_ptspec_GLRT.annotated$LS_filter <- "PASS"

    combined_ptspec_GLRT.annotated[combined_ptspec_GLRT.annotated$detected_using_size == FALSE & combined_ptspec_GLRT.annotated$DP < 20000, "LS_filter"] <- "FAIL"
    combined_ptspec_GLRT.annotated[combined_ptspec_GLRT.annotated$detected_using_size == FALSE & combined_ptspec_GLRT.annotated$DP < 20000, "ctDNA_plotting"] <- 1e-8

    combined_ptspec_GLRT.annotated$lollipop <- combined_ptspec_GLRT.annotated$LS_filter
    combined_ptspec_GLRT.annotated[combined_ptspec_GLRT.annotated$detected_using_size == TRUE & combined_ptspec_GLRT.annotated$not_detectable_dPCR == TRUE, "lollipop"] <- "non_dPCR"

    ggplot(combined_ptspec_GLRT.annotated, aes( x = reorder(SLX_barcode, ctDNA_plotting),
                                       y = (7.2 + ctDNA_plotting.log)))+
      geom_bar(stat =  "identity", width = 0.5, colour = "white")+
      scale_fill_manual(values = "chocolate1")+
      geom_point(aes(shape = lollipop), size = 3)+
      scale_shape_manual(name = "",
                         limits = c("FAIL", "non_dPCR"),
                         labels = c("low sensitivity", "non_dPCR"),
                         values = c(1, 19))+
      scale_y_continuous(breaks = c(-0.8, 0.2, 1.2, 2.2, 3.2, 4.2, 5.2, 6.2),
                    labels = c("ND", "ND", 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1))+
      theme_classic()+
      labs(x = "Sample",
           y = "IMAF",
          title = paste0(study, ": All IMAFs"))+
      theme(axis.text.x=element_blank())
    ggsave(paste0(plot_dir, "p15_waterfall_IMAF.pdf"), width = 10, height = 5)

    ## replot waterfall with fractions of cancer genomes
    # assuming the following: 1ng = 300AC or 150GE

    mutation_list_output <- read.csv(paste0(output_dir, "tumour_mutation_per_pt.csv"), header = T)
    ## get pt id
    waterfall.cancer_genomes <- left_join(combined_ptspec_GLRT.annotated, mutation_list_output[,c("Patient", "mutations")], by = "Patient")
    waterfall.cancer_genomes$cancer_genomes_fraction <- as.numeric(as.character(waterfall.cancer_genomes$mut_sum))/waterfall.cancer_genomes$mutations

    ## if non-detected, zero cancer genomes
    waterfall.cancer_genomes[waterfall.cancer_genomes$detected_using_size == FALSE,"cancer_genomes_fraction"] <- 2

    ggplot(waterfall.cancer_genomes, aes( x = reorder(SLX_barcode, ctDNA_plotting),
                                       y = cancer_genomes_fraction))+
      geom_bar(stat =  "identity", width = 0.7, colour = "white")+
      scale_fill_manual(values = "chocolate1")+
      scale_y_log10(breaks = c(1e-4,1e-3,1e-2,1e-1,1,2,10,100),
                    labels = c(1e-4,1e-3,1e-2,1e-1,1,"ND",10,100))+
      theme_classic()+
      labs(x = "Samples, ordered by IMAF",
           y = "Cancer genomes in sample")+
      theme(axis.text.x=element_blank())
    ggsave(paste0(plot_dir, "p16_waterfall_cancer_genomes.pdf"), width = 10, height = 5)

    ########## Supplementary 8a
    # 1ng = 300AC or 150GE
    # we are assuming 50% library efficiency
    # theoretical sensitivity with one locus assay would be 1/AC_input

    ## whats the proportion of detected samples that would not be detected with dPCR - 64%
    undetected <- length(which(filter(combined_ptspec_GLRT.annotated, detected_using_size)$not_detectable_dPCR))/nrow(filter(combined_ptspec_GLRT.annotated, detected_using_size))
    cutoff_line <- data.frame(c(1:80))
    colnames(cutoff_line) <- "ng_input"
    cutoff_line$cutoff <- 3.3/(cutoff_line$ng_input*300)

    filter(combined_ptspec_GLRT.annotated, adjusted_IMAF > 0, detected_using_size == "TRUE")%>%
      ggplot(aes(x = input_into_library_ng/300, y = adjusted_IMAF))+
      geom_line(data = cutoff_line, aes(x = ng_input, y = cutoff))+
      geom_point(aes(color = not_detectable_dPCR))+ #adjusted_IMAF < 3.3/(actual_lib_input*300)
      theme_classic()+
      labs(x = "Library input (ng)",
           y = "IMAF",
          title = paste0(study, ": ctDNA level of detected samples"),
          subtitle = paste0("line = 3.3/AC input; ", undetected *100, "% of samples undetected"))+
      scale_y_log10(breaks=c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1))+
      xlim(0,80)+
      scale_color_manual(breaks = c("FALSE", "TRUE"),
                         labels = c("detectable", "not detectable"),
                         name = "Detection with dPCR",
                         values = c("dodgerblue", "brown2"))
    ggsave(paste0(plot_dir, "p17_dPCR_comparison.pdf"), width = 8, height = 5)

    write.csv(waterfall.cancer_genomes, paste0(output_dir, "waterfall.cancer_genomes.csv"), row.names = F)

    ########## extract the which mutations are actually in the mut_sum . Author: Hui Zhao
    rds_files <- list.files(input_dir, pattern = ".combined.os_0.rds", full.names = T)

    rds_combined.prefilter <- data.frame()
    for (k in 1:length(rds_files)) {
      curr <- readRDS(rds_files[k])
      curr <- filter(curr, data == "ptspec")

      rds_combined.prefilter <- rbind(rds_combined.prefilter, curr)
    }

    rds_combined <- filter(rds_combined.prefilter, LOCUS_NOISE.PASS == TRUE,
                           PASS == TRUE,
                           BOTH_STRANDS == TRUE)
    ############### plotting 3bp context, using combined classes.
    context_plot <- left_join(rds_combined, SLX_layout[,c("SLX_barcode", "Timepoint")], by = "SLX_barcode")

    write.csv(context_plot[which(context_plot[,"ALT_R"]>0),],file.path(output_dir,"Mutations.for.mut_sum.csv"),row.names=F)
}

# Launch it.

if (system2('hostname', '-s', stdout = TRUE) == 'nm168s011789' && rstudioapi::isAvailable()) {
    # Rich's machine
    setwd('/home/data/INVAR')
    invisible(main(richTestOptions()))
} else {
    invisible(main(parseOptions()))
}
