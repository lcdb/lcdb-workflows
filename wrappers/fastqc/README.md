# Wrapper for FastQC

Example:

```
rule fastqc:
    input: 'samples/{sample}.fastq'
    output: 'samples/{sample}.fastqc.html'
    params:
    wrapper:
        "0.2.0/bio/fastqc"
```
