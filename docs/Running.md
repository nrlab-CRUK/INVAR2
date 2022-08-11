# Running the INVAR2 Pipeline

Assuming one has got an analysis directory set up as documented
in [Setting Up the INVAR2 Pipeline](SettingUp.md), the pipeline is ready to go.

These examples assume the _nextflow_ executable is on your path. You can check
this by typing:

    which nextflow

If you get a path back, you're ready to go. If not, you'll need to add the
directory you downloaded _nextflow_ to to your `PATH`.

## Execution Profiles

_This is covered in full [in the Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles),_
_but for this guide we'll cover only what's in the pipeline as delivered._

The INVAR2 pipeline has three "profiles" for execution called "_standard_",
"_bigserver_" and "_slurm_". These allow different types of machine or cluster
to be used without changing the pipeline code at all.

"_standard_" is the default profile that is used when no specific one is
given on the command line. It's set up for a reasonably sizeable desktop or
laptop PC. It allows a maximum of 20GB of memory and 6 cores to be used by
the pipeline. Make sure your PC will cope with this if you're using this profile.

"_bigserver_" is a profile for a powerful server with a lot of RAM and quite
a few cores. It allows a maximum of 180GB of memory and 28 cores to be used.

Both these profiles run with Nextflow's "local" executor, which will run the
pipeline's tasks on the same machine as _nextflow_. It runs as many tasks at
once as will fit into the maximum memory and number of cores specified in the
profile.

"_slurm_" is a profile for running on a Slurm cluster. It has no limit on the
maximum memory or cores; instead jobs are submitted through the Slurm job
management system asking for node space according each task. (It does limit
the maximum number of cores per job to 16, though that is by design of this
pipeline, not a limit imposed by Slurm or by Nextflow.)

__Tip__ If none of these profiles fit, you can define your own in your project's
`nextflow.config` file.
[See the Nextflow documentation.](https://www.nextflow.io/docs/latest/config.html#config-profiles)
For example:
```
profiles {
    mylaptop {
        params.MAX_CORES = 4
        process.executor = 'local'
        executor.$local.cpus = params.MAX_CORES
        executor.$local.memory = '16g'
    }
}
```

__Tip__ If you need to run on a cluster framework other than Slurm, you will
also have to define a profile for it. Each framework has its own configuration
options: [see the section on Executors in the Nextflow documentation.](https://www.nextflow.io/docs/latest/executor.html)
You should add the line `params.MAX_CORES = N` to the profile so that this
pipeline knows not to ask for more than `N` cores per job. It must be defined
for every profile.

## Starting the Pipeline

Running the INVAR2 pipeline involves just one command, run from the project's
directory:

```
nextflow run nrlab-CRUK/INVAR2
```

That is it! Nextflow will download the pipeline's files from GitHub and it's
container from DockerHub and should just work.

If you wish to use a profile other than "standard", you need to add the `-profile`
option to the command line:

```
nextflow run nrlab-CRUK/INVAR2 -profile slurm
```

_If you're running on a cluster, the command to start Nextflow should probably_
_also be submitted as a job. That is the user's responsibility._

If your project's configuration is not called `nextflow.config`, use the `-c`
option to tell Nextflow what is is called (here it's called "alternate.config"):

```
nextflow run nrlab-CRUK/INVAR2 -profile slurm -c alternate.config
```

When the pipeline starts it will check your configuration to make sure all is
well, then start processing. How long the pipeline will take depends on the
size of server it is running on or how busy the cluster it's using is.

## Changing Tasks' Resources

We've tuned the memory requirements for the processes based on the
datasets run at CRUK-CI and looking at the Nextflow execution report for
memory use. It might be that other datasets will not run with the requirements
given and will need changing. You have three main ways of altering the resources
for tasks.

### Project specific tweaks

_This is covered in full [in the Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#process-selectors)._

You can add a `process` block to your project's `nextflow.config` file to
explicitly change the memory allocated to a task or family of tasks. For example:

```Groovy
process {
    withName: createSequenceDictionary {
        memory '4gb'
    }
}
```

If you find a task that doesn't have enough memory in a particular project,
this is the way to go.

### Take a local copy for temporary use

You can clone the INVAR2 GitHub repository locally on your system, make your
changes and run from there. If you do this, the `nextflow run` command should
be given the path to the `invar.nf` file in the local clone. For example,
if your clone the repository into `$HOME/INVAR2`, the command becomes:

```
nextflow run $HOME/INVAR2/invar.nf
```

### Fork the repository for your own permanent changes

If you would like to permanently keep the changes you make to our INVAR pipeline,
you should fork the repository on GitHub into your personal GitHub area (or your
organisation's). You can then clone that fork into your working area and make the
changes you want to, pushing them back to the fork to keep them permanently and
to share with colleagues. You would run Nextflow with the name of your fork rather
than the master project.

```
nextflow run <mygithub>/INVAR2
```

If you have changes or fixes that would benefit the original code, by all means
create a pull request for us to review and incorporate where appropriate.

Forking and pull requests are beyond the scope of this document. Refer to the
numerous Git guides on the web for help.
