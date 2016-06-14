__author__ = "Justin Fear"
__copyright__ = "Copyright 2016, Justin Fear"
__email__ = "justin.m.fear@gmail.com"
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

# Figure out the prefix based on the input index, which has the format
#
#   prefix.N.bt2
#
# where N is [1-4]. We strip off the .N.bt2 and ensure the remaining prefixes
# are the same.
#
prefixes = list(set(map(lambda x: '.'.join(x.split('.')[:-2]), snakemake.output)))
assert len(prefixes) == 1, 'Multiple prefixes detected from "{0}"'.format(snakemake.output)
prefix = prefixes[0]

shell(
    "bowtie2-build "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.input} "
    "{prefix} "
    "{log}"
)
