__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
from snakemake.shell import shell


prefix = os.path.splitext(snakemake.output[0])[0]

shell("samtools view -b -q 20 {snakemake.input[0]} > {snakemake.output[0]}")
