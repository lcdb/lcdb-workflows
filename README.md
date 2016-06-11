# lcdb-workflows
snakemake workflows for LCDB bioinformatics

## First-time setup
### Conda
If you don't already have `conda`, install
[miniconda](http://conda.pydata.org/miniconda.html) (Python 3 version
recommended).

### Generate references and indexes
- Change to the `workflows/references` directory
- Create the references environment:

```bash
conda create -n references --file references_requirements.txt -c bioconda
```

- Edit `references_config.yaml` as needed. To run the tests on a different
  machine, you probably only need to edit the `data_dir` field.

- Run `workflows/references.snakefile`. This will take some time (an hour?)
  depending on connection speeds and how many CPUs are available; the minimal
  command is:

```bash
snakemake -s references.snakefile --configfile references_config.yaml
```

## Testing
There are several methods of testing. In order for a test to be run, we need
the environment to be set up, we need references to be generated (see above)
and we need example data. The test data and Snakefile are in the `test` dir;
check that directory for output.

### Method 0: development
Usually during development you'll be doing iterative changes on a local machine
or an interactive node. In this case, it's easiest to create a conda
environment based on `test/requirements.txt`:

```bash
conda create -n test-env --file test/requirements.txt -c bioconda
source activate test-env
```

Then from within the `test` directory run the example workflow using this
minimal command (you'll probably want to set `-j` as appropriate):

```bash
snakemake --configfile config.yaml
```

### Method 1: local server or interactive node
To ensure reproducibility, the most complete way to run the test is to build
a brand new environment:

```bash
# Note the trailing ".", pointing to the top-level of the repo to test
test/run_test.py --build-env .
```

This will build an environment called `lcdb-workflows-$USER-test`, deleting it
first if it already exists. It uses your current `conda` installation and
installs packages from `test/requirements.txt`.

Then run the test with:

```bash
# Note the trailing "."
test/run_test.py .
```
The test will download data (via `test/get-data.sh`) if needed and run the
`test/Snakefile` with the config file `test/config.yaml`.

If you're on an interactive node, the test will run with as many threads as
have been allocated (`$SLURM_CPUS_PER_TASK`). Otherwise, specify `--threads` in
the call to `run_test.py`.


### MEthod 2: run on cluster
Generate a submission script and print the `sbatch` command which can be pasted
to submit the job:

```bash
test/run_test.py --sbatch --threads=8 --mem=32g
```

This will also automatically generate a new environment; paste the command into
the terminal to submit the job.
