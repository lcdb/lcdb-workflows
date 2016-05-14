__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

import os
from snakemake.shell import shell
outdir = os.path.dirname(snakemake.output[0])
if not outdir:
    outdir = '.'
try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

if snakemake.log:
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""
shell(
    'multiqc '
    '--quiet '
    '--outdir {outdir} '
    '--force '
    '--filename {snakemake.output} '
    '{extra} '
    '{snakemake.params.analysis_directory} '
    '{log}'
)
