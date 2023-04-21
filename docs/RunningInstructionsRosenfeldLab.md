# Running instructions for members of the Rosenfeld lab, CRUK CI

## Set-up

Build your mutation list file and layout file as instructed in the documentation. Save both as a .csv and transfer them to the cluster ```scp -r  *.csv clust1-headnode.cri.camres.org://scratcha/nrlab/ditter01/invar_PARADIGM/``` (change this to your directory, obviously)

## Copy nextflow.config and change paths

Get a copy of the nextflow.config file that defines all your variables and file paths. There is a copy stored in /scratcha/nrlab/resources/INVAR2/ so do ```cp /scratcha/nrlab/resources/INVAR2/nextflow.config ./``` while in your directory on the cluster. 

Change the paths to where you have stored your bam files (Note they MUST be on scratcha/scratchb for the compute nodes to see them), as well as any other constants or flags you wish to use. (See [Parameters](docs/Parameters.md) for all the details.)

## Software dependencies
Make sure you have nextflow and java installed on your path:
```spack load /3tj4uay``` and ```spack load openjdk@17```

## Run it :)

In your directory with all the files, run ```nextflow run nrlab-CRUK/INVAR2 -profile slurm```

If any of your input .csv files are incorrectly filled in then red warnings will pop up and the process will stop. Else a list of processes will appear and will update as the code progresses. For a sample size of 36 patients, it takes approximately 60mins to run.




