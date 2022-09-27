# Results and Analysis

The INVAR2 pipeline produces raw outputs and a set of plots with a report from
its analysis. The outputs can be used for further analysis.

## Results Files

The raw output results of the pipeline are written to a "`results`" directory (by default,
can be changed). These are all R _RDS_ files that can be loaded back into R for
further work or conversion to other tabular formats.

The directory will contain these files.

1. `mutation_table.rds` - The on target reads that overlap the genomic regions defined in the input mutation list csv, annotated with flags indicating
all filter tests including the INVAR filter or not (the `OUTLIER.PASS` column).
No rows are filtered out from this table after selecting on target mutations: all
subsequent operations keep all rows, adding columns to increase the information.
2. `locus_error_rates.on_target.rds` - The locus error rates for the mutations
flagged as on target.
3. `error_rates.off_target.cosmic.rds` - A list of four tables giving the error
rates for off target mutations including those flagged with COSMIC.
4. `error_rates.off_target.no_cosmic.rds` - A list of four tables giving the error
rates for off target mutations excluding those flagged with COSMIC.
5. `size_characterisation.rds` - The size characterisation profile.
6. `invar_scores.rds` - The final INVAR scores.

## Plots and Analysis

The pipeline produces a series of plots and csv files written to the "`analysis`" directory (by
default). This contains an HTML report, `<STUDY>_invar2_analysis.html`, that explains the results and plots. A CSV file tracking the number of mutations per sample across the pipline is included in "`mutations_tracking.csv`", and a results summary is presented in tabular form in "`Results_summary.csv`".
