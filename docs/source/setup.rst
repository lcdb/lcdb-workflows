First-time setup
================

1. Install `conda`.

`lcdb-workflows` depends heavily on the `conda` package manager to ensure
reproducible environments. If you don't already have `conda`, install
`miniconda <http://conda.pydata.org/miniconda.html>`_ (Python 3 version
recommended).

2. Clone the repository::

   git clone https://github.com/lcdb/lcdb-workflows.git

3. Run the tests. Set `--threads` to as many cores as you have on your machine.
With 8 cores, this takes about 6 min::

   test/run_test.py workflows/mapping/Snakefile --clean --build-env --threads 8

The output will be in `workflows/mapping/pasilla`.

Subsequent setup
================

The `lcdb/lcdb-init.py` script can be used to initialize a working repo in
a new working directory.
