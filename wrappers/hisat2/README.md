# Wrapper for HISAT2

[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) performs graph-based
alignment of next generation sequencing reads to a population of genomes.


## Input

* `fastq`: either a single FASTQ file (for single-end) or a list of two FASTQ
  files (for paired-end). HISAT2 will auto-detect compression based on
  extenstion (`.gz` or `.bz2`)

* `index`: List of the `*.ht2` index files. See the example below for a good
  way to create this list. The wrapper will figure out the prefix to provide to
  HISAT2 based on these files.


## Output

* Sorted BAM file of aligned reads

## Threads

Threads are passed to HISAT2 as the `--threads` parameter.

## Params

* `extra` can be a string of arbitrary parameters that will be passed verbatim
  to HISAT2.

## Example

This example runs HISAT2 on single-end reads, reducing the max intron length
and modifying the output for use with downstream transcriptome assemblers:

```python
rule hisat2:
    input:
        fastq="mapped/{sample}_R1.fastq.gz",
        index=expand('/data/assembly/assembly.{n}.ht2', n=range(1,9))
    output:
        "mapped/{sample}.sorted.bam.bai"
    params: extra="--max-intronlen=50000 --dta"
    wrapper:
        "/path/to/hisat2/wrapper.py"
```
