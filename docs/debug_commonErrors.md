# A list of Common errors and their solutions

### 1 - Column `UNIQUE_POS` doesn't exist, followed by "There are no rows for sample XXX in mutation_table.on_target.rds".

You probably have different patient names in the layout file and in the mutation list (check those O/0's) so when the on target loci are chosen the PATIENT_CHROM:POS columns don't match 
