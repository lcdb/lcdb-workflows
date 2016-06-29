Testing
=======

The interface for testing is through the `test/run_test.py` script. Its
built-in help should give you all you need to know:

.. code:: bash

    test/run_test.py -h


Here are some common use-cases. These examples use the "mapping" workflow.

**Minimal test:** Assumes environment `lcdb-workflows-$USER-env` exists and data
have been downloaded. If you're on an interactive node, all your allocated
cores will be used. This is the version you'll probably use the most when
devleoping a new wrapper or rule or workflow.

.. code:: bash

    test/run_test.py workflows/mapping/Snakefile


**Only rebuild environment:** Made some changes to requirements.txt and want to
refresh your env? Build a new one without running the test:

.. code:: bash

    test/run_test.py workflows/mapping/Snakefile --build-env --no-test


**Build environment, clean output,  and run test:** A more complete test, this
destroys an existing environment and builds it again. It then **deletes the
entire test data directory** and re-downloads test data, then runs the test:

.. code:: bash

    test/run_test.py workflows/mapping/Snakefile --build-env --clean

**Submit job to a single cluster node:** This runs a complete test on a single cluster node.
It builds a new environment before submitting to the cluster. Use the
`--threads` and `--mem` args to configure the job. An email will be sent to you
upon job completion or failure.

.. code:: bash

    test/run_test.py workflows/mapping/Snakefile --sbatch --build-env --clean --threads=8 --mem=32g

**Submit job in parallel to cluster:** This is probably what you will be doing
in "production" mode -- running jobs in parallel, each submitted to a different
node. Note that the config file contains per-rule thread and memory settings,
so `--threads` and `--mem` are unused.

.. code:: bash

    test/run_test.py workflows/mapping/Snakefile--cluster --clean --build-env
