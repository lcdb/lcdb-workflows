# Wrapper for picard MarkDuplicates.

[Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
identifies duplicate reads.

## Input

Coordinate sorted and indexed bam file.

## Output

A new bam file with duplicate reads marked in the SAM flags field, or
optionally removed.

## Log

Log output will be saved to specified log file, otherwise it will go to stdout

## Threads

Ignores threads.

## Params

Additional parameters can be passed to picard verbatim by supplying a string in
params.extra.

# Example

```
rule picard_MarkDuplicates:
    input:
        bam = "mapped/{sample}.sorted.bam",
        bai = "mapped/{sample}.sorted.bam.bai"
    output:
        bam = "mapped/{sample}.sorted.dup.bam"
        metrics = "mapped/{sample}.sorted.dup.metrics"
    params:
        extra="REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT"
    wrapper:
        "/path/to/picard/MarkDuplicates/wrapper.py"
```
