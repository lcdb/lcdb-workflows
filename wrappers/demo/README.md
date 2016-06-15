# Wrapper for *addtoolnamehere*


## Input
Describe expected inputs. Be sure to point out if they are expected to have
some property (e.g., if an input BAM needs to be coord sorted, name sorted, or
indexed; if fastqs are required to be paired-end, etc).

If the wrapper expects inputs to be named, specify that here. Here's an
example from, say, an aligner:

* `fastq`: fastq file; optional compression detected by extension (.gz or .bz2)
* `index`: paths to bowtie2 `*.bt2` index files. The index prefix will be
  inferred from these paths.

## Output
Describe output file[s].

## Threads
Mention if multiple processors can be utilized. Provide the optimal number of
processors to use if known.

## Params
* Describe any parameters to be provided via the `params:` keyword and how they
  are used by the wrapped program.
* Typically at least `extra` should be provided, which can be passed verbatim
  to the programs.
* Try to give some pointers on what to use.

## Example
* Include the example Snakemake rule in triple-backtick format, specifying
  python as the language (to get nice syntax highlighting on github).

* As a convention, provide the wrapper path as `/path/to/wrapper/program`, and
  include the `{sample}` placeholder in input and output (if applicable).

* Try to provide more than just a minimal example, since it's easier to delete
  stuff than to add.

```python
rule samtools_sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    params: extra="-n"
    wrapper:
        "/path/to/wrapper/location"
```
