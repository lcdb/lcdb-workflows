__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

try:
    extra = snakemake.params['extra']
except AttributeError:
    extra = ""

if snakemake.log:
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""

# Handle paired-end reads. Since snakemake automatically converts a one-element
# list to a string, here we detect single-end reads by checking if input.fastq
# is a string.
if isinstance(snakemake.input.fastq, str):
    fastqs = '-U {0} '.format(snakemake.input.fastq)
else:
    assert len(snakemake.input.fastq) == 2
    fastqs = '-1 {0} -2 {1} '.format(*snakemake.input.fastq)

# Figure out the prefix based on the input index, which has the format
#
#   prefix.N.ht2
#
# where N is [1-8]. We strip off the .N.ht2 and ensure the remaining prefixes
# are the same.
#
prefixes = list(set(map(lambda x: '.'.join(x.split('.')[:-2]), snakemake.input.index)))
assert len(prefixes) == 1, 'Multiple prefixes detected from "{0}"'.format(snakemake.input.index)
prefix = prefixes[0]

shell(
    "hisat2 "
    "-x {prefix} "
    "{fastqs} "
    "--threads {snakemake.threads} "
    "{extra} "
    "-S {snakemake.output}.sam "
    "{log}"
)

shell("samtools view -Sb {snakemake.output}.sam > {snakemake.output} && rm {snakemake.output}.sam")
