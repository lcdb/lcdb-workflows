# Wrapper for samtools view

## Input

SAM or BAM.

## Output

SAM or BAM, depending on parameters

## Params

* `extra` is passed verbatim to `samtools view`


## Example:

Remove multimappers from a BAM file, creating a new BAM file:

```
rule samtools_sort:
    input: "mapped/{sample}.bam"
    output: "mapped/{sample}.unique.bam"
    params: extra = "-b -q 20"
    wrapper:
        "samtools/view"
```
