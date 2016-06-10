__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell
try:
    extra = params['extra']
except KeyError:
    extra = ""

if len(snakemake.input.fastq) == 1:
    fastqs = '-U {snakemake.input}'
else:
    assert len(snakemake.input.fastq) == 2
    fastqs = '-1 {snakemake.input.fastq[0]} -2 {snakemake.input.fastq[1]}'

shell(
    "hisat2 "
    "-x {params.index} "
    '{fastqs} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.output}.sam"
)

shell("samtools view -Sb {snakemake.output}.sam > {snakemake.output}")

