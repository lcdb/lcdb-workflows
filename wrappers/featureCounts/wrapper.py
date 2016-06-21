__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

# Supports arbitrary arguments to be passed in to the shell call without
# raising an error if no params.extra were provided.
try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

# Redirects stdout and stderr to a log if it is provided.
if snakemake.log:
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""

# creates counts.full and counts.full.summary
shell(
    "featureCounts -T {snakemake.threads} "
    "{extra} "
    "-a {snakemake.input.annoGTF} "
    "-o {snakemake.output.countsFull} "
    "{snakemake.input.mappingFile} "
    "{log} ")

# subsets counts.full to include meta-feature w/counts
shell(
    "tail -n +3 {snakemake.output.countsFull} "
    "| cut -f 1,7 > {snakemake.output.counts} "
)
