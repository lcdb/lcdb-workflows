__author__ = "Your Name"
__copyright__ = "Copyright 2016, Your Name"
__email__ = "your@email.edu"
__license__ = "MIT"

# include any imports that will be needed, example below
from snakemake.shell import shell

# include commands with params, input and output from snakemake rule, compare README.md with the example below
shell("samtools index {snakemake.params} {snakemake.input[0]} {snakemake.output[0]}")
