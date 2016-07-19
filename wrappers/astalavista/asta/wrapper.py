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
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""

outfile = snakemake.input.gtf + "_astalavista.gtf.gz"
shell(
    "astalavista -t asta "
    "-i {snakemake.input.gtf} "
    "{extra} "
    "{log} ")
shell("mv {outfile} {snakemake.output}")
