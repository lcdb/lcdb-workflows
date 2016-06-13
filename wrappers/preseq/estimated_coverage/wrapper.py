__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

# import snakemake's ability to execute shell commands
from snakemake.shell import shell

# execute bam2mr
shell("bam2mr {snakemake.input[0]} > sorted.nr")

# execute preseq gc_extrap
shell("preseq gc_extrap -o {snakemake.output[0]} sorted.nr")

# remove sorted.nr
shell("rm sorted.nr")
