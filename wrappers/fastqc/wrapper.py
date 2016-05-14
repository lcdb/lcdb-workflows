__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

import os
from snakemake.shell import shell
outdir = os.path.dirname(snakemake.output[0])

# fastqc creates a zip file and an html file but the filename is hard-coded by
# replacing fastq|fastq.gz|fq|fq.gz|bam with _fastqc.zip|_fastqc.html in the
# input file's basename.
#
# So we identify that file and move it to the expected output after fastqc is
# done.

outfile = os.path.basename(snakemake.input[0])
strip = ['.fastq', '.fq', '.gz', '.bam']
for s in strip:
    outfile = outfile.replace(s, '')
out_zip = os.path.join(outdir, outfile + '_fastqc.zip')
out_html = os.path.join(outdir, outfile + '_fastqc.html')

shell(
    'fastqc '
    '--threads {snakemake.threads} '
    '--noextract '
    '--quiet '
    '--outdir {outdir} '
    '{snakemake.params.extra} '
    '{snakemake.input} '
)
shell('mv {out_zip} {snakemake.output.zip}')
shell('mv {out_html} {snakemake.output.html}')
