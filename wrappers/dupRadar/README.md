# Wrapper for dupRadar

dupRadar provides an easy way to distinguish between artifacts vs natural
duplicate reads in RNA-Seq data. 

[Link to homepage](https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html)

[Link to manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)

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
rule dupRadar:
    input:
       dupBAM = '{sample_dir}/{sample}/{sample}.cutadapt.hisat2.unique.sort.dedup.bam',
       GTF = '/data/LCDB/references/dm6/dm6_r6-11.gtf.gz',
       dupRScript = wrapper_for('/dupRadar/dupRadar.R')
    output:
        dupSmooth = '{sample_dir}/{sample}/{sample}_dupRadar_drescatter.png'
    params:
        paired="no",
        stranded="no"
    wrapper:
        wrapper_for('dupRadar')
```
