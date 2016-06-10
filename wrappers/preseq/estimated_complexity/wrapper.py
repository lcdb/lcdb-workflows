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
    
# execute preseq c_curve
shell("preseq lc_curve {snakemake.params.BAM} {snakemake.input[0]} {snakemake.params.output} {snakemake.output[0]}")
