#!/bin/bash

# Runs specified workflows, building their rulegraph and full DAG, saving them
# as PNGs in the docs image directory.

set -eo pipefail

source activate lcdb-workflows-$USER-env

OUTDIR=source/images
mkdir -p $OUTDIR

# Optional args for dot formatting
RULEGRAPH_DOT=""

# arrange graph left-to-right since these can get pretty unwieldy
DAG_DOT="-Grankdir=LR"
workflow=$1
config=$2
type=$3

if [ $type == "rulegraph" ]; then
    snakemake \
        --rulegraph \
        --forceall \
        --configfile $config \
        --directory ../workflows/$workflow \
        -s ../workflows/$workflow/Snakefile \
    | dot -Tpng $RULEGRAPH_DOT > $OUTDIR/${workflow}_dag.png
fi

if [ $type == "dag" ]; then 
    snakemake \
        --dag \
        --forceall \
        --configfile $config \
        --directory ../workflows/$workflow \
        -s ../workflows/$workflow/Snakefile \
    | dot -Tpng $DAG_DOT > $OUTDIR/${workflow}_full_dag.png
fi
