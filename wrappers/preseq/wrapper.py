__author__= "Behram Radmanesh"
__copyright__= "Copyright 2016, Behram Radmanesh"
__email__= "radmaneshbs@nih.gov"
__license__= "NIH"

"""Work In Progress"""

import os
from snakemake.shell import shell
outdir = os.path.dirname(snakemake.output[0])

# Handle BAM files
if (snakemake.input[0] == {sample}.bam):
    shell(
        samtools sort {sample}.bam -o {sample}.sorted.bam
        preseq c_curve -o complexity_output.txt -B {sample}.sorted.bam
        preseq lc_extrap -o future_yield.txt -B {sample}.sorted.bam
        bam2mr {sample}.sorted.bam > sorted.bam.nr
        gc_extrap -o future_coverage.txt sorted.bam.nr
    )
else (snakemake.input[0] == {sample}.bed):
    shell(
        sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 {sample}.bed > {sample}.sort.bed
        preseq c_curve -o complexity_output.txt {sample}.sort.bed
        preseq lc_extrap -o future_yield.txt {sample}.sort.bed
        samtools sort {sample}.bam -o {sample}.sorted.bam
        bam2mr {sample}.sorted.bam > {sample}.sorted.nr
        preseq gc_extrap -o future_coverage.txt {sample}.sorted.nr
    )
