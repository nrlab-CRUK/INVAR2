# Results and Analysis

The INVAR2 pipeline produces raw outputs and a set of plots with a report from
its analysis. The outputs can be used for further analysis.

## Results Files

The results of the pipeline are written to a "`results`" directory (by default,
can be changed). These are all R _RDS_ files that can be loaded back into R for
further work or conversion to other tabular formats.

The directory will contain these files.

1. `mutation_table.rds` - The on target mutations, annotated with flags indicating
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

The pipeline will write the analysis results into an "`analysis`" directory (by
default). This will contains an HTML report, `<STUDY>_invar2_analysis.html`,
and numerous plots in PDF format. It will also contain some summary CSV files from
its own processing.
