Testing
=======
Provide `run_test.py` a directory pointing to the top level of a repo. Data
will be downloaded if needed into that directory and the test will run in that
directory (so make sure you have enough room).

For now the tests use data from test/get-data.sh, but once we have test data
in place they will run on that. The tests also use the already-prepared
references/indexes/annotations in /data/LCDB/references.

Use it like this when actively developing on an interactive node. It will use
as many cores as you have allocated (using `$SLURM_CPUS_PER_TASK`), and will
use whatever environment you are currently in:

```bash
run_test.py /path/to/repo
```

As a shortcut to building the necessary environment, run it like this to build
a conda environment called lcdb-workflows-$USER-test using the requirements
(deleting the env if it already exists). You'll need to source it afterwards:

```bash
run_test.py --build-env /path/to/repo
```

"The works": force --build-env, and print out the command to submit to the
cluster using the specified --threads and --mem.

```bash
run_test.py --sbatch /path/to/repo
```
