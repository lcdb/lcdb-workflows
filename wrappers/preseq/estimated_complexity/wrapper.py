__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

from snakemake.shell import shell
try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

shell("preseq lc_extrap {extra} -B {snakemake.input[0]} -o {snakemake.output[0]}")
