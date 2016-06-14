# Wrapper for Bowtie2 align

Bowtie2 is a basic alignment tool for short reads. This wrapper is for running
bowtie2 alignments.

[Link to homepage](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[Link to manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

## Input

* `fastq`: either a single FASTQ file (for single-end) or a list of two FASTQ
  files (for paired-end). HISAT2 will auto-detect compression based on
  extenstion (`.gz` or `.bz2`)

* `index`: List of the `*.bt2` index files. See the example below for a good
  way to create this list. The wrapper will figure out the prefix to provide to
  Bowtie2 based on these files.

*Note:* If pair-end reads are used, be sure to add the correct mate pair
directs to `params.extra`, in other words add `--fr/--rf/--ff`. For Illumina
pair-end reads us `--fr`.

## Output

* Sorted BAM file of aligned reads

## Threads

Threads are passed to Bowtie2 as the `--threads` parameter.

## Params

* `extra` can be a string of arbitrary parameters that will be passed verbatim
  to Bowtie2.

# Example

This example runs Bowtie2 on single-end reads:

```
rule bowtie2:
    input: fastq='samples/{sample}.fastq',
           index=expand('/data/assembly/assembly.{n}.bt2', n=range(1,4))
    output: 
            "mapped/{sample}.bt2.sorted.bam.bai"
    params: extra="--local --very-fast-local"
    log: samples/{sample}.bt2.log
    wrapper:
        "/path/to/bowtie2/align/wrapper.py"
```
