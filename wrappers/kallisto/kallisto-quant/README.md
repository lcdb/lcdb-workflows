# Wrapper for *kallisto quant*


[Kallisto home](https://pachterlab.github.io/kallisto/)

[Kallisto manual](https://pachterlab.github.io/kallisto/manual)


## Input

* `index`: kallisto index file
* `fastq`: either a pair of PE fastqs or a single fastq. If single, see note in Params section.

## Output

`kallisto quant` outputs files hard-coded as:

```
outfile/
├── abundance.h5
├── abundance.tsv
└── run_info.json
```

Since the downstream `sleuth` tool expects these names, we do not move them. It
is recommended that you use the sample name in the directory (see example
below). The wrapper detects the output directory from the first output file;
the rest of the output files are optional.

## Threads

Threads are passed as the `--threads` argument to kallisto quant.

## Params

* `extra`: passed verbatim to kallisto quant

NOTE: if you are using single-end reads, then kallisto will complain if you
don't set the `--fragment-length` and `--sd` params. Add them to extra. Based
on [this
comment](https://groups.google.com/forum/#!searchin/kallisto-sleuth-users/sd/kallisto-sleuth-users/VPJfzL502bw/e2JDq7ezBgAJ)
from the author, `--fragment-length=200 --sd=20` should be a good starting
point.

NOTE: if you plan on using `sleuth` for differential expression, you'll want to
include the `--bootstrap-samples` argument.

## Example

Single-end reads example, with all output files specified (only one is required
for the wrapper to identify the output directory), log, and `extra` params:

```python
rule kallisto_quant:
    input:
        fastq="{sample}.fastq.gz",
        index="references/dm6/kallisto.idx"
    output:
        "quant/kallisto/{sample}/abundance.h5",
        "quant/kallisto/{sample}/abundance.tsv",
        "quant/kallisto/{sample}/run_info.json"
    log: "quant/kallisto/{sample}/kallisto.log"
    params: extra="--fragment-length=200 --sd=20 --single"
    wrapper:
        "/path/to/wrapper/location"
```

Paired-end, simplified example:

```python
rule kallisto_quant:
    input:
        fastq=["{sample}_R1.fastq.gz", "{sample}_R2.fastq.gz"]
        index="references/dm6/kallisto.idx"
    output:
        "quant/kallisto/{sample}/abundance.h5",
    wrapper:
        "/path/to/wrapper/location"
```
