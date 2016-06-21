# Wrapper for dupRadar

dupRadar provides an easy way to distinguish between artifactual vs natural
duplicate reads in RNA-Seq data. Prior to dupRadar only global duplication rates
were used and they don't take into account the effect of gene expression levels. 
dupRadar relates *duplication rates* and *length normalized read counts* of every
gene to model the dependency of both variables. 

[Link to homepage](https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html)

[Link to manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)

## Input
* `BAM`: BAM file with mapped reads has to be duplicate marked using either Picard
  or BamUtil
* `GTF`: GTF file contaning the genetic features to count the reads falling on the
  features.

## Output
* `PNG file`: an expression density scatter plot in .png format is generated

## Threads
Threads not supported.

## Params
* `paired`: for paired end sequencing, default is *no*
* `stranded`: for stranded reads, default is *no*

## Example

```python
rule dupRadar:
    input:
       dupBAM = '{sample_dir}/{sample}/{sample}.cutadapt.hisat2.unique.sort.dedup.bam',
       GTF = '/data/LCDB/references/dm6/dm6_r6-11.gtf.gz',
       dupRScript = wrapper_for('/dupRadar/dupRadar.R')
    output:
        dupSmooth = '{sample_dir}/{sample}/{sample}_dupRadar_drescatter.png'
    params:
        paired="no",
        stranded="no"
    wrapper:
        wrapper_for('dupRadar')
```
