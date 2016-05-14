# Wrapper for cutadapt

[cutadapt] trims adapters and poly-A tails and has flexible options for quality
and hard-trimming.

## Input
FASTQ

## Output
Trimmed FASTQ file, and optionally other files depending on the provided params.

## Log
Log output will be saved to specified log file, otherwise it will go to stdout

## Threads
Ignores threads

## Params
Additional parameters can be passed to cutadapt verbatim by supplying a string in params.extra.

# Example

```
rule cutadapt:
    input: 'samples/{sample}.fastq'
    output: 'samples/{sample}.trim.fastq'
    params: extra="-a file:adapters.fa -q 20"
    wrapper:
        "file://path/to/cutadapt"
```
