#!/bin/sh

if [[ -d spython ]]
then
    source spython/bin/activate
else
    python3 -m venv spython
    source spython/bin/activate
    pip install spython
fi

spython recipe Dockerfile > singularity_spec.txt

sudo rm -rf invar_sandbox
sudo singularity build --sandbox invar_sandbox singularity_spec.txt


# See https://stackoverflow.com/a/60316979
# sudo singularity build my_container.sif docker-daemon://local/my_container
