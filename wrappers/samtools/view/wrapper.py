__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

if snakemake.log:
    log = "2> {log}".format(snakemake.log)
else:
    log = ""

shell("samtools view {extra} {snakemake.input[0]} > {snakemake.output[0]} {log}")
