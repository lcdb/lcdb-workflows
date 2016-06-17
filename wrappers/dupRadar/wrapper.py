__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

if snakemake.log:
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""

shell(
    "Rscript {snakemake.input.dupRScript} "
    "{snakemake.input.dupBAM} "
    "gtf={snakemake.input.GTF} "
    "stranded={snakemake.params.stranded} "
    "paired={snakemake.params.paired} "
    "outfile={snakemake.output[0]} "
    "threads={snakemake.threads} {snakemake.params.name}"
)
