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

1. The reference genome FASTA file, with index.
3. [1,000 genomes](https://www.internationalgenome.org/) SNP VCF file.
4. [COSMIC](https://cancer.sanger.ac.uk/cosmic/) coding mutants VCF file.

It also needs two files that need to be curated for the studies that use INVAR2: the layout
file and the tumour mutations file for the study being processed.

### Layout File

This is a CSV file containing information about all studies and the samples in them. The
columns the INVAR2 pipeline requires are listed below. Other columns may be present but
play no part in the INVAR2 processing. Neither is column order important.

    STUDY
    SAMPLE_ID
    BAM_FILE
    CASE_OR_CONTROL
    PATIENT
    INPUT_INTO_LIBRARY_NG
    SAMPLE_NAME
    SAMPLE_TYPE
    TIMEPOINT

`PATIENT` is the patient id and is used to link to rows in the tumour mutations file.
`SAMPLE_ID` can be any text to internally identify a sample: it just must be unique
in this file.
`CASE_OR_CONTROL` identifies cases (cancer patients) or control samples, and its content
should either be "`case`" or contain "`control`", such as "`control_negative`". This
file can contain information for all projects, so the `STUDY` column provides the means of
listing samples for multiple data sets.

`BAM_FILE` is the name of the aligned BAM file for this sample. It will be relative
to the `BAM_PATH` parameter, which by default will be a directory called `bam` in your
project directory. One can add an optional column `ACTIVE` to the layout file to prevent
some BAM files being used for the analysis. If this column is present and has a value that
is not either "true" or "yes", then the sample will not be included in the analysis. If the
column is present but there is no value for a sample, it's interpreted as being included.
In other words, if you explicitly want a sample excluded, put "no" or "false" into the column
for the sample.

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

Note: please ensure each loci is called only once per patient! (ie no duplicates in the list) else this would cause issues in teh sizeAnnotation.R step.

## Creating a Project Directory

You should create a directory structure for your analysis. This forms the top level
directory under which you should place your data files (aligned BAM) and the other files
we will need to create for INVAR2 to run. As the pipeline runs, other files and directories
will be created containing the results.

__Tip__ This is the directory where you will run _nextflow_ from. In the Nextflow files, this
directory is accessed by the built in variable "`launchDir`", e.g. `"${launchDir}/bam"`.

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

### BAM Files

The simplest way to make your BAM files available to the INVAR2 pipeline is to
copy them into a directory `bam` in the project directory. This though can mean
making copies of large files that could be accessed in their original locations,
if those locations are visible to the machine(s) you are running INVAR2 on.
In these circumstances, one can set the
parameter "`BAM_PATH`" in your project's `nextflow.config` file
to list all the directories your BAM files can be found in ([see the parameters page](Parameters.md)).
The pipeline will then be able to read the files from their original locations without
making copies.

__Tip__ Be wary of using symbolic links. Much of the INVAR2 pipeline runs inside Singularity and
links that point outside the regular file structure may not be visible to the pipeline's
tasks. This can cause errors that might surprise, as outside of Singularity those files
will be accessible.

__Tip__ If a file cannot be found, the pipeline will stop with an error that's a
little odd. For example, if your file is called "SampleA.bam" but the pipeline
can't find it in any directory given by `BAM_PATH`, you will see an error such as:

    No such file: /dev/null/SampleA.bam

Ignore the "_/dev/null_" bit: this is just a way to force the pipeline to stop
when a BAM file cannot be found. What it's really saying is that the file "SampleA.bam"
cannot be found in any directory listed in `BAM_PATH`.
