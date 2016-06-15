__author__ = "Jusitn fear"
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

shell(
    "picard MarkDuplicates "
    "I={snakemake.input.bam} "
    "O={snakemake.output.bam} "
    "{extra} "
    "M={snakemake.output.metrics} "
    "{log}"
)
