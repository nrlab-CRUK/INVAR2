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
    HG19_GENOME = "${REFERENCE_ROOT}/hg19.genome"
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
| STUDY                     |      | The name of the study being processed.                     |
| ERROR_SUPPRESSION_NAME    |      | Setting giving the description of the error suppression mechanism used in processing the BAM files. e.g. "f0.9_s2" |
| FAMILY_SIZE               | int  | The family size.                                           |
| LAYOUT_TABLE              |      | The path to the layout table file.                         |
| TUMOUR_MUTATIONS_CSV      |      | The path to the patient tumour mutations file.             |
| HG19_GENOME               |      | Path to the hg19 genome index file.                        |
| FASTA_REFERENCE           |      | Path to the hg19 FASTQ reference file.                     |
| THOUSAND_GENOMES_DATABASE |      | Path to the 1,000 genomes SNP file (VCF).                  |
| COSMIC_DATABASE           |      | Path to the COSMIC variants file (VCF).                    |

## Optional Parameters

These parameters all have defaults defined in the `nextflow.config` file in the
pipeline's source directory. Any redefined in your project's `nextflow.config`
file will take precedence over the defaults.
Types are strings unless indicated by `int`eger, decimal `num`ber or `bool`ean (`true`/`false`).

| Parameter                      | Type | Default                   | Description/Purpose                                       |
|--------------------------------|------|---------------------------|-----------------------------------------------------------|
| INPUT_FILES                    |      | "${launchDir}/to_run.csv" | Path to the `to_run.csv` file listing source BAM files to include in the analysis. If you wish to name this file something else, or have it in a different location, you can change this parameter. |
| MAPPING_QUALITY                | int  | 40                        | Minimum mapping quality threshold.                        |
| BASE_QUALITY                   | int  | 20                        | Minimum base quality threshold.                           |
| MPILEUP_MINIMUM_DEPTH          | int  | 2                         | Minimumin depth to consider for mpileup. Set to 1 for sWGS samples. |
| MPILEUP_MAXIMUM_DEPTH          | int  | 100000                    | Maximum depth for mpileup.                                |
| SLOP_BASES                     | int  | 10                        | How many bases either side of the target base to assess for the background error rate. |
| REMOVE_DUPLICATES              | bool | `true`                    | Whether to remove duplicates in pile ups.                 |
| MAXIMUM_DEPTH                  | int  | 1500                      | Omit data points with uncharacteristially high unique depth given the input mass used. |
| MINIMUM_REFERENCE_DEPTH        | int  | 5                         | Here we require at least 5 reference reads at a locus. Set to 0 for sWGS. |
| MQSB_THRESHOLD                 | num  | 0.01                      | Exclude data points due to poor MQ and SB.                |
| ALT_ALLELES_THRESHOLD          | int  | 3                         | Blacklist loci with &ge; N separate alternate alleles.    |
| MINOR_ALT_ALLELE_THRESHOLD     | int  | 2                         | Blacklist multiallelic loci with a mutant read count of &ge; N in the minor mutant allele. |
| COSMIC_THRESHOLD               | int  | 0                         | Loci with &gt; N entries in COSMIC are considered as COSMIC mutations. |
| PROPORTION_OF_CONTROLS         | num  | 0.1                       | Blacklist loci that have signal in &gt; P of the non-patient specific samples. |
| MAXIMUM_BACKGROUND_MEAN_ALLELE_FREQUENCY | num | 0.01             | Filter loci with a background allele frequency in controls greater than this value. |
| ALLELE_FREQUENCY_THRESHOLD     | num  | 0.01                      | Maximum allele frequency value for acceptable samples.    |
| IS_BLOODSPOT                   | bool | `false`                   | Only change to true if you are running blood spot data through the pipeline. This omits outlier-suppression on samples with deduplicated depth of &lt;5x because high AF loci cannot be reliably identified with low depth. |
