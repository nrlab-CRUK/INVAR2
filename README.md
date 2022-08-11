# INVAR2
#Paulius Mennea PhD student is working on this branch to enable INVAR2 on sWGS data

A Nextflow based software tool to run the INVAR (Integration of VAriant Reads)
analysis pipeline. Previously available [here.](https://bitbucket.org/nrlab/invar/wiki/Home)

The pipeline seeks to detect minimal residual disease signal in patient liquid biopsy data,
and outputs a classifcation based on a threshold score established on a cohort of healthy
samples. A general cohort size should range from around 20+ individual patients. 

More details are available in the links below. Test data is available at the EGA accession numbers EGAS00001004447 and EGAS00001005246
with results available in the following publications:

[1] Wan J, Heider K, Gale D et al. ctDNA monitoring using patient-specific sequencing and integration of variant reads. Sci Transl Med. 2020;12(548). doi:10.1126/scitranslmed.aaz8084

[2] Heider K, Wan J, Gale D. ctDNA detection by personalised assays in early-stage NSCLC. MedRxiv. doi:https://doi.org/10.1101/2021.06.01.21258171


Contact: Hui.zhao@cruk.cam.ac.uk, emma-jane.ditter@cruk.cam.ac.uk

## Documentation

1. [Setting up the INVAR2 pipeline.](docs/SettingUp.md)
2. [INVAR2 parameters.](docs/Parameters.md)
3. [Running INVAR2.](docs/Running.md)
4. [Results and Analysis.](docs/ResultsAndAnalysis.md)
5. [Indepth Technical Details.](docs/TechnicalDetails.md)
6. [Running Instructions for Rosenfeld Lab members.](docs/RunningInstructionsRosenfeldLab.md)
7. [FAQ's](docs/TechnicalNotes_FAQ.md)
