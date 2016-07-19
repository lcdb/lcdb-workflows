#!/bin/bash
set -e
if [ -n "$SLURM_JOBID" ]; then
    export TMPDIR="/lscratch/$SLURM_JOBID"
fi
