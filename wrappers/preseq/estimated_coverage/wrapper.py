__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

import os
# import snakemake's ability to execute shell commands
from snakemake.shell import shell

# specify an output directory from snakemake otherwise insert default
outdir = os.path.dirname(snakemake.output[0])
if not outdir:
    outdir = '.'

try:
    output = snakemake.params.output
except AttributeError:
    output = ""


# execute bam2mr
shell("bam2mr {snakemake.input[0]} > sorted.nr")

# execute preseq gc_extrap
shell("preseq gc_extrap {snakemake.params.quick} {snakemake.params.paired} -o {snakemake.output[0]} sorted.nr")

# remove sorted.nr
#shell("rm sorted.nr")
