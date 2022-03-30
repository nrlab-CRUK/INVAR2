# In depth details on INVAR2 functioning


The INVAR2 pipeline is run through invar.nf, available to the top level of the GitHub directory.

Before INVAR2 begins processing the input data, high level software checks are made to ensure the input files can be located, the layout table is in the correct format, and all parameters have been defined.

All the genomic locations that are in the patient panels are extracted from the input bam files and are annotated with COSMIC SNV and 1000 genome SNP information. The remaining loci then pass through 3 filters, LOCUS_NOISE.PASS that removes noisy loci with too many alt alleles or mutated in too many samples, BOTH_STRANDS.PASS, which makes sure the alt allele is found on both the forward and reverse strand, and OUTLIER_SUPPRESSION.PASS, which removes loci who's signal is an outlier. The remaining loci are then used to estimate the quantity of ctDNA in a sample, which is used to calculate the log likelihood of a mutation being present in the sample. The lo likelihood value is used to calculate an INVAR score, and detection is determined based on a threshold of 95\% specificity on the healthy control samples. 

More details in the working of each script is given below. 

## Step 1

The process begins by in "INVAR2/processes/1_parse.nf"

The bed file (tumour mutation file) is processed with slop in "INVAR2/templates/1_parse/slopPatientInfo.sh", and the bam files are processed such that only the regions defined in the bed file (tumour mutation file) is conserved. The pileups are done in "INVAR2/templates/1_parse/mpileup.sh"
The multiallelic sites are split into biallelic records and the columns of interest are conserved ("INVAR2/templates/1_parse/biallelic.sh").
A dataframe containing every base pair in the pileup (here +/- 10bp from the loci defined in the tumour mutations file) is created and annotated with the tabix information defining if it is a SNP as known from the tabix database (tabix dataframe used). ("INVAR2/templates/1_parse/tabix.sh").
Another dataframe with COSMIC mutation information is created for the same genomic positions.

Returning to the raw data, a dataframe containing every read in the supplied bam files that overlaps with the genomic regions defined in the bed file/tumour mutations file is then annotated with all this information and more (tabix, cosmic, trinucleotide) in the "INVAR2/python/1_parse/addTabixAndTrinucleotides.py" script.

The dataframe above is then filtered by the defined MQSB threshold and blacklists loci that have more than 3 alt alleles or insufficient/too much depth. It also annotates dataframe with True/False values if the loci passes the COSMIC and 1000 genomes thresholds. Done by "INVAR2/R/1-parse/createMutationsTable.R", it outputs a hidden mutation_table.filtered.rds file (to find the location type ```find -name mutation_table.filtered.rds``` while in the INVAR2 directory).

The first of three filter flags (LOCUS_NOISE.PASS) is defined in "INVAR2/R/1-parse/offTargetErrorRates.R" where the error rates of each loci are then calculated in "INVAR2/R/1-parse/offTargetErrorRates.R". This mean the information of each base pair on each read is condensed into one row in the dataframe that corresponds to a single genomic location in the tumour mutations file (ie the read information is condensed into a single set of metrics). The loci error rates dataframes are saved by COSMIC on_target or off_target flag in "INVAR2/Results/error_rates.off_target.*.rds".


Calls mutationTable which is 






## Step 2

## Step 3

## Step 4

## Step 5

The script "GeneralisedLikelihoodRatioTest.R" calculated the loglikelihood of ctDNAbeing present in each sample, and outputs it as an INVAR score. It calculated an INVAR score for when the fragment size distribution is taken into account, and an INVAR score for when the fragment sizes are neglected. The script reads in INVAR2/results/mutation_table.rds and calculates the distribution of fragment lengths of all samples and removes the fragments belonging to the sample in question (to avoid a circular calculation). It calculates the probability distribution of reads for healthy and mutant reads in the range of the predefined read lengths ( minFragmentLength to maxFragmentLength as defined in the parameters). There is the posibility of weighting non-mutant molecules to a constant. \textcolor{red}{WHY IS IT 0.1?}. The probabilities of each length are kept as probabilities if the read is mutant, and set to 0.1 if not. \textcolor{red}{WHY ARE THERE MUTANT READS IN HEALTHY DISTRIBUTION?!!! line 205 in generaliseLikelihoodRatioTest.R excludes mutant. Does it mean to exclude "cases"?}.




