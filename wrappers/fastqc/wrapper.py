__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

import os
from snakemake.shell import shell
outdir = os.path.dirname(snakemake.output[0])
shell(
    'fastqc ',
    '--threads {snakemake.threads} '
    '--noextract '
    '--outdir {outdir} '
    '{snakemake.input} '
)
