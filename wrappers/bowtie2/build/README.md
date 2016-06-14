# Wrapper for Bowtie2 Build

# Wrapper for Bowtie2 build

Bowtie2 is a basic alignment tool for short reads. This wrapper is for running
bowtie2-build to create bowtie2 indices.

[Link to homepage](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[Link to manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

## Input

A reference file in FASTA format.

## Output

Generates a set of index files for use with bowtie2 aligner.

## Threads

Supports threads, passed in as the `--threads` arg

## Params

Additional parameters can be passed to verbatim by supplying a string in `params.extra`.

# Example

```
rule bowtie2build:
    input: "{prefix}.fasta"
    output: expand("bowtie2/{{prefix}}.{N}.bt2", n=[1,2,3,4])
    params: extra="--seed 123"
    threads: 2
    log: "bowtie2/{prefix}.bt2.build.log
    wrapper:
        "/path/to/bowtie2/build/wrapper.py"
```
