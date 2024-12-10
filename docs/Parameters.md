# INVAR2 Parameters

Here are all the parameters one can set to control the INVAR2 pipeline.
They need to be set in your project's `nextflow.config` file in a `params`
block. For example:

```
params {
    STUDY = "MyStudy"
    ERROR_SUPPRESSION_NAME = "f0.9_s2"
    FAMILY_SIZE = 2
    ...
}
```

This file is (almost) a [Groovy](http://groovy-lang.org/) file, and the
language's rules apply. String values must be enclosed in single or
double quotes. Numbers should be written without any quotes. Logical
(boolean) values must use the Java literals `true` or `false` (no quotes).
No parameter in the INVAR2 pipeline needs lists or maps so fortunately
we don't need to worry about those.

__Tip__ There is a difference in Groovy between strings enclosed in single
and double quotes. Single quotes are literal strings as they are written.
Double quotes are what are called _GStrings_ in the Groovy world, and can
refer to other parameters using the `${}` syntax. For example, if the
reference data is in the same directory, one could write:

```
params {
    REFERENCE_ROOT = "/home/reference_data/INVAR2"
    FASTA_REFERENCE = "${REFERENCE_ROOT}/ucsc.hg19.fasta"
    COSMIC_DATABASE = "${REFERENCE_ROOT}/CosmicCodingMuts.vcf.gz"
}
```

It is safe to add parameters for your own convenience such as `REFERENCE_ROOT`
in this example. The pipeline won't know anything about them directly but won't
object to their presence.

## Required Parameters

These parameters must be defined in the project specific `nextflow.config`
file as there is no possible default for them. Parameters are all upper case.
Types are strings unless indicated by `int`eger, decimal `num`ber or `bool`ean (`true`/`false`).

| Parameter                 | Type | Description/Purpose                                        |
|---------------------------|------|------------------------------------------------------------|
| STUDY                     |string| The name of the study being processed.                     |
| ERROR_SUPPRESSION_NAME    |string| Setting giving the description of the error suppression mechanism used in processing the BAM files. e.g. "f0.9_s2". The variable is used in the title of the output plots, allowing users to easily compare two INVAR runs with differently processed bam files. |
| FAMILY_SIZE               | int  | The family size used when deduplicating the BAM file. Used to annotate output plots and is not used in numerical calculation of results.                                           |
| LAYOUT_TABLE              |      | The path to the layout table file.                         |
| TUMOUR_MUTATIONS_CSV      |      | The path to the patient tumour mutations file.             |
| FASTA_REFERENCE           |      | Path to the reference genome FASTQ file.                   |
| THOUSAND_GENOMES_DATABASE |      | Path to the 1,000 genomes SNP file (VCF).                  |
| COSMIC_DATABASE           |      | Path to the COSMIC variants file (VCF).                    |

## Optional Parameters

These parameters all have defaults defined in the `nextflow.config` file in the
pipeline's source directory. Any redefined in your project's `nextflow.config`
file will take precedence over the defaults.
Types are strings unless indicated by `int`eger, decimal `num`ber or `bool`ean (`true`/`false`).
The default settings are for targeted sequencing. For Whole Genome Sequencing, you may want to adjust the parameters. 
Some parameters are dependent on the sequencing depth (DP) of your plasma bam files. Recommanded settings are listed.

| Parameter                      | Type | Default                   | Description/Purpose                                       |
|--------------------------------|------|---------------------------|-----------------------------------------------------------|
| BAM_PATH                       |      | "${launchDir}/bam"        | A list of directories the aligned BAM files can be found in. |
| RESULTS_DIR                    |      | "${launchDir}/results"    | The directory to write results (R _RDS_ files) to.        |
| ANALYSIS_DIR                   |      | "${launchDir}/analysis"   | The directory to write analysis plots and report to.      |
| MAPPING_QUALITY                | int  | 40                        | Minimum mapping quality threshold.                        |
| BASE_QUALITY                   | int  | 20                        | Minimum base quality threshold.                           |
| MPILEUP_MINIMUM_DEPTH          | int  | 2                         | Minimum depth to consider for mpileup. Set to 1 for sWGS samples. |
| SLOP_BASES                     | int  | 10                        | How many bases either side of the target base to assess for the background error rate. |
| REMOVE_DUPLICATES              | bool | true                      | Whether to remove duplicates in pile ups.                 |
| MAXIMUM_DEPTH                  | int  | 1500                      | Maximum depth to remove abnormally high depth loci. Recommand to set to 10*DP. |
| MINIMUM_REFERENCE_DEPTH        | int  | 5                         | Minimum number of reference reads of a locus to be considered. Recommand to set to 1. Set to 0 for sWGS. |
| MQSB_THRESHOLD                 | num  | 0.01                      | Exclude data points due to poor MQ and SB.                |
| ALT_ALLELES_THRESHOLD          | int  | 3                         | Blacklist loci with &ge; N separate alternate alleles.    |
| MINOR_ALT_ALLELE_THRESHOLD     | int  | 2                         | Blacklist multiallelic loci with a mutant read count of &ge; N in the minor mutant allele. |
| COSMIC_THRESHOLD               | int  | 0                         | Loci with &gt; N entries in COSMIC are considered as COSMIC mutations. |
| PROPORTION_OF_CONTROLS         | num  | 0.1                       | LOCUS_NOISE.PASS filter: If the mutations in the patient are also found in &gt; P number of other samples, then the mutations would be filtered out at this step. If only using case samples, recommand to set to at least `2 / total cases`. If there are control samples, then set to `1~2 / total controls`. |
| MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY | num | 0.01             | Filter loci with a background allele frequency in controls greater than this value. |
| ALLELE_FREQUENCY_THRESHOLD     | num  | 0.01                      | Maximum allele frequency value for acceptable samples. This value affects the CONTAMINATION_RISK.PASS filter by preventing samples with too high AF be control samples. Recommand to set to 0.1.   |
| MAXIMUM_MUTANT_READS           | int  | 10                        | Maximum number of mutant reads of a locus to be considered for detection. Recommand to set to 1*DP. |
| MINIMUM_INFORMATIVE_READS      | int  | 20000                     | Minimum number of informative reads (reads that have passed all filters) for each sample to be considered as part of the cohort. Ie a low sensitivity threshold (samples with fewer than MINIMUM_INFORMATIVE_READS after filtering) will not be considered for classification. Recommand to set to 100*DP.                      |
| IS_BLOODSPOT                   | bool | false                     | Only change to true if you are running blood spot data through the pipeline. This omits outlier suppression on samples with deduplicated depth of &lt;5x because high AF loci cannot be reliably identified with low depth. |
| OUTLIER_SUPPRESSION_THRESHOLD  | num  | 0.05                      | Outlier suppression threshold. Related to                  |
| MINIMUM_FRAGMENT_LENGTH        | int  | 60                        | Minimum fragment length.                                  |
| MAXIMUM_FRAGMENT_LENGTH        | int  | 300                       | Maximum fragment length.                                  |
| SMOOTHING                      | num  | 0.25                      | Smoothing function for size profile (width of smoothing). |
| ONLY_WEIGH_MUTANTS             | bool | true                      | Only weigh ctDNA signal based on mutant fragments.        |
| SCORE_SPECIFICITY              | num  | 0.95                      | Score specificity for ROC plot.                           |
| ITERATIONS                     | int  | 10                        | Number of iterations when subsampling the wild-type reads (control samples) to determine a range of INVAR scores. |
| LAYOUT_TABLE_ENCODING          |      | "UTF-8"                   | Character encoding for the layout file.                   |
| TUMOUR_MUTATIONS_CSV_ENCODING  |      | "ASCII"                   | Character encoding for the mutations file.                |

`BAM_PATH` is a list of directories in which the BAM files can be found. It is
handled very much like the _`PATH`_ environment variable one is familiar with
if using Unix. It is a list of directory paths separated by a colon (semicolon in
Windows). Each directory is examined in turn to try
to find each BAM file listed in the layout file; the first directory in which the
file appears is used for the full path to the file for the pipeline.
Thus one can specify more than one directory to look for BAM files in.
The pipeline will halt with an error if a file cannot be
found in any of the directories given by `BAM_PATH`.
