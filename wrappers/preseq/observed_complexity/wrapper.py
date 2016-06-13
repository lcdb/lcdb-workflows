__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

# import snakemake's ability to execute shell commands
from snakemake.shell import shell

# execute preseq c_curve
shell("preseq c_curve -B {snakemake.input[0]} -o {snakemake.output[0]}")
