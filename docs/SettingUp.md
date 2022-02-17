# Setting Up the INVAR2 Pipeline

## Prerequisites

The INVAR2 pipeline runs within the [Nextflow pipeline system](https://www.nextflow.io/).
Instructions for installing Nextflow on your system can be found in
[the Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html).

You will also need Java installed on your system. Nextflow is written in
[Groovy](http://groovy-lang.org/), which is a scripting language written on top of Java.
There is no Groovy installation, but the Java system is required as a separate install.

You should put the `nextflow` executable somewhere on your `PATH`. If Java isn't installed
in a standard location, you will also need to add the `java` executable to your `PATH` too.

## Reference Data

The INVAR2 pipeline needs reference data to perform its processing. Some of this is publicly
available standard reference data, and some is study specific information created before
running INVAR2.

The pipeline needs these publicly available files:

1. An index file for the UCSC hg19 human genome.
2. The UCSC hg19 reference FASTA file, with index.
3. [1,000 genomes](https://www.internationalgenome.org/) SNP VCF file.
4. [COSMIC](https://cancer.sanger.ac.uk/cosmic/) coding mutants VCF file.

It also needs two files that need to be curated for the studies that use INVAR2: the layout
file and the tumour mutations file for the study being processed.

### Layout File

This is a CSV file containing information about all studies and the samples in them. The
columns the INVAR2 pipeline requires are listed below. Other columns may be present but
play no part in the INVAR2 processing. Neither is column order important.

    STUDY
    SAMPLE_NAME
    SAMPLE_TYPE
    CASE_OR_CONTROL
    POOL
    BARCODE
    INPUT_INTO_LIBRARY_NG
    PATIENT
    TIMEPOINT

`PATIENT` is the patient id and is used to link to rows in the tumour mutations file.
`POOL` and `BARCODE` are ids for the pool a sample is in and the barcode used in sequencing
for that sample respectively; they cross link into the `to_run.csv` file discussed below.
`CASE_OR_CONTROL` identifies cases (cancer patients) or control samples, and its content
should either be "`case`" or contain "`control`", such as "`control_negative`". This
file can contain information for all projects, so the `STUDY` column provides the means of
listing samples for multiple data sets.

### Tumour Mutations File

This is a CSV file listing patient mutations. The columns it must have are listed below.
Again, order is unimportant and other columns may exist in the file.

    CHROM
    POS
    REF
    ALT
    TUMOUR_AF
    PATIENT

`CHROM` is the chromosome [e.g. "chr1"], `POS` the location on the chromosome of the
mutation, `REF` the reference allele, `ALT` the alternate allele, `TUMOUR_AF` the
tumour allele frequency [a decimal], `PATIENT` the patient identifier. `PATIENT` is
cross linked to the `PATIENT` column in the layout file.


## Creating a Project Directory

You should create a directory structure for your analysis. This forms the top level
directory under which you should place your data files (aligned BAM) and the other files
we will need to create for INVAR2 to run. As the pipeline runs, other files and directories
will be created containing the results.

__Tip__ This directory where you will run _nextflow_ from. In the Nextflow files, this
directory is accessed by the built in variable "`launchDir`", e.g. `"${launchDir}/bam"`.

### BAM Files

For simplicity, copy the BAM files into a directory `bam` in the top level directory.
If these files are already available on the file system, one could hard link them rather
than copy them to save space. If they are on another area of the file system, they needn't
be copied but their full path needs to be listed in the `to_run.csv` file discussed below.

__Tip__ Be wary of using symoblic links: much of the INVAR2 pipeline runs inside Singularity and
links that point outside the regular file structure may not be visible to the pipeline's
tasks. This can cause errors that might surprise, as outside of Singularity those files
will be accessible.

### `to_run.csv`

This file is a summary of the BAM files. It needs to list the pool, barcode and file name
for each BAM that should be part of the analysis. Files can be removed from this list
to exclude them from processing without having to remove them from the `bam` directory.

The file must start with a header line that contains exactly the values (any order):

    POOL
    BARCODE
    FILE_NAME

The `POOL` and `BARCODE` columns cross linked to the same columns in the layout file.
We would expect then each row to form a unique pair in this file and a unique pair within
a study in the layout file. If a pool + barcode pair does not exist in the layout file,
the pipeline will stop with an error.

`FILE_NAME` gives the name of the BAM file for the pool + barcode. It can be a relative
path to the top level directory (e.g. "`bam/sample1.bam`") or an absolute path.

### Configuration File

The configuration file sets the parameters for the INVAR2 run. There are some parameters
that must be defined in it, and others that can be set to override the defaults in the
INVAR2 pipeline as provided.

This file should be called "`nextflow.config`" and its structure is described in the
[Nextflow documentation](https://www.nextflow.io/docs/latest/config.html). One can put
settings in this file for everything described in that web page, but for typical use one
need only define a `params` block.

__Tip__ The file doesn't have to be called `nextflow.config`, but if it is it will be
picked up automatically when the pipeline is run. One could have more than one configuration
file in the directory and choose which one to use with _nextflow_'s `-c` option on the
command line.

Typically the file will look like:

```
params {
    <parameter1> = <value>
    
    <parameter2> = <value>
}
```

It will have as many parameters as needed. [The parameters page](Parameters.md) lists
all the parameters that can or must be set to run INVAR2.
