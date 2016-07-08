__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

import os
from snakemake.shell import shell

try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

if snakemake.log:
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""

outdir = os.path.dirname(snakemake.output[0])

shell(
    "kallisto quant "
    "--index {snakemake.input.index} "
    "-o {snakemake.output} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.input.fastq} "
    "{log} ")
