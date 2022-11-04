# In depth details on INVAR2 functioning

A stream-of-though document to store Frequently Asked Questions and important technical information which might be useful for the wider community.

## Running FAQ's

###  "CHROM column missing in mutation list" error

Make sure you export your mutation list to "Comma-seperated Values (.csv)" NOT "CSV UTF-8 format"

### Do you have a recommended 1000 Genome SNP file?
We have used the 1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz file. Some users have experienced issues with templates/1_parse/tabix.sh (error message: Failed to read BGZF block data at offset) switching to the recommended file has sorted it. 


## Technical FAQ's

### What is ON_TARGET variable?

Defined in createMutationsTable.R, line 88: ON_TARGET = UNIQUE_POS %in% tumourMutationsTable$UNIQUE_POS
so whether it is a loci as defined in the tumour mutation file or not (ie if the mutation is at the desired loci or at a different location in the read)?

### How is the final classification done?

* Invar scores are given for each sample as called against each panel, and are the log likelihoods of there being mutations in the sample given the signal at the pre-chosen loci
* Ideally, we would measure very low scores for non patient specific measurements (ex: patient A's data as called against patient B's panel, ie low levels of background noise for randomly chosen loci) and high scores when measuring signal across the patient specific panel (ie a high prevalence of the tumour mutations detectable in the plasma)
* A threshold value is selected to give our results a specificity (true negative) score of eg. 95%. 
* The threshold score is used to classify each sample based on the sample's patient specific INVAR score. 

### In Results_summary.csv, what is MUTATION_SUM and MUTATIONS?
MUTATION_SUM = sum(MUTANT) , so number of mutated reads but with no filtering conditions.  

MUTATIONS = n_distinct(PATIENT, UNIQUE_POS) so number of input mutations or loci of interest.

