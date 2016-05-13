# lcdb-workflows
snakemake workflows for LCDB bioinformatics

[![Build Status](https://travis-ci.org/lcdb/lcdb-workflows.svg?branch=master)](https://travis-ci.org/lcdb/lcdb-workflows)

## Testing
Continuous integration testing is performed by
[Travis-CI](https://travis-ci.org) on every push to github.

See the `.travis.yml` file for the configuration. The testing takes advantage
of a prepared [Docker](https://docker.com) container that has all prerequisites
installed, and then runs the `travis-test.sh` script on the travis-ci infrastructure.

The Docker container has already been created and uploaded to
[dockerhub](https://hub.docker.com) using the commands:

```bash
cd docker && docker build -t daler/smklo . && docker upload daler/smklo
```


### Run local tests with Docker
This ensures the most isolated environment to run local tests, and most closely
approximates the testing performed on travis-ci. If the tests pass locally,
they should pass on travis-ci. But first you need to do some setup:

1. Install Docker ([Mac OSX](https://docs.docker.com/mac/) or
   [Linux](https://docs.docker.com/linux/)).

2. Download the Docker container with:

```bash
docker pull daler/smklo
```

3. Download and prepare the example data:

```bash
test/get-data.sh
```

This completes the setup. In the future, to run the tests:

```bash
docker run --rm -it -v $(pwd):/opt/lcdb -u $(id -u):$(id -g) daler/smklo /bin/bash travis-test.sh
```

Some explanation of what's happening there:

- `docker run` runs the container. Think a virtual machine that starts up instantly
- `--rm` removes the container when it exits to save space
- `-it` means `--interactive` `--tty`. It attaches the command line to the container after starting it so it essentially runs in the foreground
- `-v $(pwd):/opt/lcdb` exports the current working directory into the
  container, effectively mounting it at `/opt/lcdb`. The container can then make
  changes to the directory and they will show up on the filesystem.
- `-u $(id -u):$(id -g)` sets the user and group of the container to the
  current user and group. This is to avoid having the container create files
  owned by root that you'd have to later chown.
- `daler/smklo` is the name of the container
- `/bin/bash travis-test.sh` runs bash as the shell and then runs the test driver script.

If you want to poke around in the container, leave off the `travis-test.sh` and
you'll drop into a shell in the container. Note that any changes you make to
the system (apt-get install, conda install, etc) will not persist after the
container exits.

## Setting up without Docker
Not all environments support Docker (i.e. biowulf/helix). In this situation,
you can use the `docker/requirements.txt` file to build your own isolated
environment.

If you don't already have `conda`, install
[miniconda](http://conda.pydata.org/miniconda.html) (Python 3 version
recommended).

Then choose a name for your environment. Here we use `lcdb`:

```bash
conda create -n lcdb --file requirements.txt -c bioconda -c r python=3
```

This creates an isolated environment. To use it you have to temporarily activate it:

```bash
source activate lcdb
```

When you're done, deactivate

```bash
source deactivate
```

Once your environment is activated, run the tests:

```bash
bash travis-test.sh
```
