# Wrapper for samtools sort.

## Example:

```
rule samtools_sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    wrapper:
        "samtools/view"
```
