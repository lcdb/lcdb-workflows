__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

# bam2mr converts BAM to a format used by preseq
shell("bam2mr {snakemake.input[0]} > sorted.nr")

try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

shell("preseq gc_extrap {extra} -o {snakemake.output[0]} sorted.nr")

# the nr is only used internally for this wrapper
shell("rm sorted.nr")
