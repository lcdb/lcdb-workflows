#!/bin/bash

# Intended to be called from the top-level dir of the repo. Installs non-conda
# dependencies, then tries to run the test snakefile.

cd test
snakemake
