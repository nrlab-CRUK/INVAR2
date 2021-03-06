# In depth details on INVAR2 functioning


The INVAR2 pipeline is run through invar.nf, available to the top level of the GitHub directory, with inputs and constants defined in the nextflow.config text-like document. 

Before INVAR2 begins processing the input data, high level software checks are made to ensure the input files can be located, the layout table is in the correct format, and all parameters have been defined.


## Pipeline Overview

All the genomic locations in the patient panels are extracted from the input bam files and are annotated with COSMIC SNV and 1000 genome SNP information. The remaining loci then pass through 3 filters: LOCUS_NOISE.PASS removes noisy loci with too many alt alleles or mutated in too many samples, and BOTH_STRANDS.PASS makes sure the alt allele is found on both the forward and reverse strand (Step 1). The data frame is annotated with the size information of all the reads (Step 2) and then the final filter OUTLIER_SUPPRESSION.PASS, or PASS, removes loci who's signal is an outlier (Step 3). The strength of the outlier suppression step is defined in the nextflow.config input "OUTLIER_SUPPRESSION_THRESHOLD", with default 0.05. The remaining loci are then used to estimate the quantity of ctDNA in a sample, which is used to calculate the log likelihood of a mutation being present in the sample. This is output as an INVAR score per sample (Step 4). The final classification of each sample (cancer detected/non detected) is determined based on a threshold of 95\% specificity on the healthy control samples (Step 5). 

More details in the working of each script is given below. 

## Step 1

This step parse's the input bam files to select the regions defined in the mutation file, create the pile up's and annotates each read with COSMIC, 1000genome and trinucleotide information. It is run in "INVAR2/processes/1_parse.nf".
It creates a data frame where each read (r_1, or forward read, and r_2, or reverse read) feature on each row, and are annotated with the trinucleotide context, COSMIC information at that location, and information from the 1000genomes project if the position is a known SNP. The error rates at each locus are output into a .rds file (results/locus_error_rates.on_target.rds) as are the off target, ie not patient specific genomic region error rates in results/error_rates.off_target.[no_/]cosmic.rds. 

The bed file (tumour mutation file) is processed with slop in "INVAR2/templates/1_parse/slopPatientInfo.sh", and the bam files are processed such that only the regions defined in the bed file (tumour mutation file) is conserved. The pileups are done in "INVAR2/templates/1_parse/mpileup.sh"
The multiallelic sites are split into biallelic records and the columns of interest are conserved ("INVAR2/templates/1_parse/biallelic.sh").
A data frame containing every base pair in the pileup (here +/- 10bp from the loci defined in the tumour mutations file) is created and annotated with the tabix information defining if it is a SNP as known from the tabix database (tabix data frame used). ("INVAR2/templates/1_parse/tabix.sh").
Another data frame with COSMIC mutation information is created for the same genomic positions.

Returning to the raw data, a tsv dataframe containing every read in the supplied bam files that overlaps with the genomic regions defined in the bed file/tumour mutations file is then annotated with all this information and more (tabix, cosmic, trinucleotide) in the "INVAR2/python/1_parse/addTabixAndTrinucleotides.py" script.

The data frame above is then filtered by the defined MQSB threshold and blacklists loci that have more than 3 alt alleles or insufficient/too much depth. It also annotates data frame with True/False values if the loci passes the COSMIC and 1000 genomes thresholds. Done by "INVAR2/R/1-parse/createMutationsTable.R", it outputs a hidden mutation_table.filtered.rds file per input sample (to find the location type ```find -name mutation_table.filtered.rds``` while in the INVAR2 directory).

The first of three filter flags (LOCUS_NOISE.PASS) is defined in "INVAR2/R/1-parse/offTargetErrorRates.R" where the error rates of each loci are then calculated in "INVAR2/R/1-parse/offTargetErrorRates.R". The information of each base pair on each read is condensed into one row in the data frame that corresponds to a single genomic location in the tumour mutations file (ie the read information is condensed into a single set of metrics). The on target locations (patient specific genomic locations) are saved in locus_error_rates.on_target.rds. For off target loci (control sample loci or non patient specific case loci) the background error rates are calculated for COSMIC flag equals True or False, and saved in "INVAR2/Results/error_rates.off_target.*.rds". 


## Step 2

This step annotates the data frame create above with the size information of each read, and summarises the size distribution of all reads into a data frame. This process is run from  "INVAR2/processes/2_size_annotation.nf". 

The step starts by creating a temporary file mutationlist.bed using the INVAR2/R/1_parse/patientListCsvToBed.R script from Step 1. This is then piped into INVAR2/templates/2_size_annotation/getFragmentSize.sh which runs INVAR2/python/2_size_annotation/getFragmentSize.py. getFragmentSize.py analyses each read in the pileup and analyses it with the alt allele and the size of the fragment as calculated by pysam version 0.17.0. Finally, the information is collated with the main data frame by INVAR2/R/2_size_annotation/sizeAnnotation.R


## Step 3

This step performs outlier suppression on all loci based on the outlier suppression constant defined in nextflow.config file. The default value is 0.05 and helps define the last of three filter: the PASS filter. 

INVAR2/R/3_outlier_suppression/outlierSuppression.R imports the osMutationsFile's from hidden directory in work/ and calculates the quantity of circulating tumour DNA in the sample through the expectation maximisation algorithm. Using the estimated quantity of ctDNA, the binomial probability of seeing N mutant reads is calculated. The final filter, the PASS filter, is then defined as the probability of seeing a mutation greater than the outlier suppression constant (default 0.05) divided by the number of loci in that sample (multiple hypothesis testing). 

## Step 4

This step performs the generalised likelihood ration test on each sample and outputs an invar score. 

The script "GeneralisedLikelihoodRatioTest.R" calculated the log likelihood of ctDNAbeing present in each sample, and outputs it as an INVAR score. It calculated an INVAR score for when the fragment size distribution is taken into account, and an INVAR score for when the fragment sizes are neglected. An INVAR score is calculated for each sample when called against the their own patient specific panel, and an invar score is calculated when other patient samples are called against the patient specific panel, to act as controls. 

The script reads in the individual mutation_table.outliersuppressed.*.rds from each sample passed through 3-outlier_suppression.nf (these are intermediate files hidden within /work/ and can be found by '''find -name mutation_table.outliersuppressed.*.rds). It then calculates the distribution of fragment lengths of all samples and removes the fragments belonging to the sample in question (to avoid a circular calculation). 

It calculates the probability distribution of reads for healthy and mutant reads across all samples (only includes case reads) in the range of the predefined read lengths (minFragmentLength to maxFragmentLength as defined in the parameters). There is the posibility of weighting non-mutant molecules to a constant. \textcolor{red}{WHY IS IT 0.1?}. The probabilities of each length are kept as probabilities if the read is mutant, and set to 0.1 if not. 

This is then fed into the calc_likelihood_ratio_with_RL function which calculated the log likelihood of ctDNA being present in the sample given an estimate of the ctDNA from the expectation maximisation algorithm. This outputs a constant, which is output as the INVAR score for that sample. Each patient panel is called against the patient sample, and a patient specific INVAR score is output. As a measure of controls, in addition to the healthy samples, each individual patient specific panel is called against other patient data, which theoretically should produce no signal. Thus, INVAR2 subsamples the reads of the other patient data and outputs an INVAR score for each subsample. This subsampling and INVAR score calculation step is repeated 10 times. 

## Step 5

This step determines the final classification of the samples and produces the plots available in analysis/ . 

The classification step is done by the function ```adjustInvarScores```line 118 of analysisCalculations.R. It uses all non patient specific INVAR scores (healthy samples and the 10 iterations of the non-patient specific data) to determine the threshold INVAR score that would give a specificity of 95\% on said data. This threshold score is then used on the patient specific INVAR scores the classify if cancer is detected in the sample or not. 
The plots in analysis/ are produced by INVAR2/R/5_analysis/analysis.R and calls INVAR2/R/5_analysis/analysisCalculations.R and INVAR2/R/5_analysis/analysisPlots.R

A .Rmd file is produced to explain the outputs in greater details. 

Step 5 is open to further personalisation depending on individual use cases. All modifications of this final step requires the whole of the pipeline to be rerun. 


Note: any feedback or suggested improvements that might be useful to the wider community are appreciated. Please email: emma-jane.ditter@cruk.cam.ac.uk, hui.Zhao@cruk.cam.ac.uk or rosenfeld@cruk.cam.ac.uk. 


