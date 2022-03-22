# In depth details on INVAR2 functioning


The INVAR2 pipeline is run through invar.nf, available to the top level of the GitHub directory.

Before INVAR2 begins processing the input data, high level software checks are made to ensure the input files can be located, the layout table is in the correct format, and all parameters have been defined.

## Step 1

The process begins by in "INVAR2/processes/1_parse.nf"

The bed file (tumour mutation file) is processed with slop in "INVAR2/templates/1_parse/slopPatientInfo.sh", and the bam files are processed such that only the regions defined in the bed file (tumour mutation file) is conserved. The pileups are done in "INVAR2/templates/1_parse/mpileup.sh"
The multiallelic sites are split into biallelic records and the columns of interest are conserved ("INVAR2/templates/1_parse/biallelic.sh").
A dataframe containing every base pair in the pileup (here +/- 10bp from the loci defined in the tumour mutations file) is created and annotated with the tabix information defining if it is a SNP as known from the tabix database (tabix dataframe used). ("INVAR2/templates/1_parse/tabix.sh").
Another dataframe with COSMIC mutation information is created for the same genomic positions.

Returning to the raw data, a dataframe containing every read in the supplied bam files that overlaps with the genomic regions defined in the bed file/tumour mutations file is then annotated with all this information and more (tabix, cosmic, trinucleotide) in the "INVAR2/python/1_parse/addTabixAndTrinucleotides.py" script.

The dataframe above is then filtered by the defined MQSB threshold and blacklists loci that have more than 3 alt alleles or insufficient/too much depth. It also annotates dataframe with True/False values if the loci passes the COSMIC and 1000 genomes thresholds. Done by "INVAR2/R/1-parse/createMutationsTable.R", it outputs a hidden mutation_table.filtered.rds file (to find the locaton type ```find -name mutation_table.filtered.rds``` while in the INVAR2 directory).

The error rates of each loci are then calculated in "INVAR2/R/1-parse/offTargetErrorRates.R" 
Calls mutationTable which is 






## Step 2

## Step 3

## Step 4

## Step 5
