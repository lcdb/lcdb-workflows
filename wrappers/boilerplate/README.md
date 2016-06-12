# Wrapper for *addtoolnamehere*

Include a brief description about the tool here.

Discuss anything about the wrapper that may be surprising or non-intuitive to
experienced users of the command-line version of the wrapped tool.

[Link to homepage](http://example.com)

[Link to manual](http://example.com)

GUIDELINES (delete this section in the actual README)

Use [GitHub-flavored
Markdown](https://help.github.com/articles/basic-writing-and-formatting-syntax/)
for formatting.

When making choices about what inputs to require, aim for simplicity in the
rules. For example, if a tool needs an indexed BAM, note that requirement in
the "Input" section but don't require both the BAM and the BAI files as input
-- just the BAM is fine.

It's too much work and frankly too brittle to wrap every command line
argument as a possible `params:` option. Instead, when possible, simply provide
an `extra` param, which can hold a string to be passed verbatim to the wrapped
tool. That way, even if program arguments change across versions the wrapper
will not have to change.

Aligners with indexes: for now, provide all indexes as inputs so that they can
be created if needed. So for bowtie2, the input indexes in the snakefile rule
should be along the lines of `expand("{prefix}.{n}.bt2", n=[1, 2, 3, 4])` and
the wrapper should expect a list.

Some tools have hard-coded output filenames. To avoid restricting
filenames over on the snakefile side, wrappers for such tools should let the
tool create its hard-coded output in a temporary directory and then move it to
the final filename[s] provided by the rule.

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
