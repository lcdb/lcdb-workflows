# Wrapper for *addtoolnamehere*

Include a brief description about the tool here.

Discuss anything about the wrapper that may be surprising or non-intuitive to
experienced users of the command-line version of the wrapped tool.

[Link to homepage](http://example.com)

[Link to manual](http://example.com)

*GUIDELINES (delete this section in the actual README)

Use [GitHub-flavored
Markdown](https://help.github.com/articles/basic-writing-and-formatting-syntax/)
for formatting the README. Line wrap at 80 characters (`set textwidth=80` in
vim). At minimum, include the sections shown below -- Input, Output, Params,
Threads, and Example.

The goal of the README should be to show how to use the wrapper and to minimize
the time it takes from discovering the wrapper to writing a rule for it. No
need to re-write or copy the manual. However if the wrapper does something
substantially different than the corresponding command line call, then spend
time to clearly document that. An example of this might be a peak-caller
warpper that outputs a bigBed file even though the peak-caller itself only
creates a BED file. Whether or not this otherwise unexpected behavior is up for
debate.

The goal of the wrapper should be to contain the complexity of running
a tool while exposing a simple interface -- input, output, and params -- to the
snakefile. Sometimes the wrapper will be barely more than a `shell()` call (see
samtools index); no need to make it any more complex. Others can get a little
involved in order to simplify how it is used (see fastq_screen).

When making choices about what inputs to require, for a wrapper, aim for
simplicity in the rules. For example, if a tool needs an indexed BAM, note that
requirement in the "Input" section of the README but don't require both the BAM
and the BAI files as input -- just the BAM is fine, since the BAI filename can
be inferred from it.

It's too much work and frankly too brittle to wrap every command line
argument as a possible `params:` option. Instead, when possible, simply provide
an `extra` param, which can hold a string to be passed verbatim to the wrapped
tool. That way, even if program arguments change across versions, the wrapper
will not have to change.

Aligners with indexes: for now, provide all indexes as inputs so that they can
be created if needed. So for bowtie2, the input indexes in the snakefile rule
should be along the lines of `expand("{prefix}.{n}.bt2", n=[1, 2, 3, 4])` and
the wrapper should expect a list. From that list, the wrapper can figure out
the prefix to provide to bowtie2. This allows the rule to specify the
dependency of the index (re-running if the index has been updated) without the
extra redundancy of also specifying an index. See the hisat2 wrapper for an
example of this.

Some tools have hard-coded output filenames. To avoid restricting
filenames over on the snakefile side, wrappers for such tools should let the
tool create its hard-coded output in a temporary directory and then move it to
the final filename[s] provided by the rule. See the fastqc wrapper for an
example of this.*

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
rule featureCounts:
    input:
		annoGTF="{sample}.gtf",
		mappingFile="{sample}.bam"
    output:
        counts="{sample}_counts.txt"
    params: 
		extra="-T 5 ",
		featureType="exon ",
		attributeType="gene_id "
    wrapper:
        "/path/to/wrapper/location"
```
