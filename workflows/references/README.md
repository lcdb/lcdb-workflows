# References workflow

The references workflow is intended to be a universal workflow that supports
arbitrary genomes. The config file allows very flexible specification of how to
create the files.

It's probably easiest to show an example config and then describe what's
happening. This example sets up a complete set of Drosophila references and
indexes, and just a bowtie2 index for E. coli (for use by fastq_screen to check
for contamination.

```yaml
dm6:
  fasta:
    url: "url to fasta file"
    postprocess: "dm6.fasta_postprocess"
    indexes:
      - bowtie2
      - hisat2
      - star
  gtf:
    url: "url to GTF.gz"
    postprocess: "dm6.gtf_postprocess"

  transcriptome:
    url: "url to transcriptome fa.gz"
    indexes:
      - kallisto

ecoli:
  fasta:
    url: "url for E. coli gfa.gz"
    indexes:
      - bowtie2
```

Top-level keys here describe assembly names (here, dm6 and ecoli). Under that
are 3 optional keys: `fasta`, `gtf`, and `transcriptome`. Each of these has at
least a `url` entry. They can also have a `postprocess`, which is an arbitrary
function (described below) that converts the downloaded URL to something that
conforms to the standards of the workflow (also described below). The `fasta`
and `transcriptome` keys can have an optional `indexes` section  which will
build the specified indexes.

## Post processing
The `postprocessing` values are dotted names that refer to Python modules
importable by the `reference.snakefile`. The dotted name should refer to
a function that has the function signature:

```python
def func(downloaded_filename, postprocessed_filename)
```

The job of a postprocessing function is to ensure that the
fastq/gtf/transcriptome fasta meets the requirements below and is ready for any
intended downstream tasks. For example if we download the fasta file from
FlyBase for dm6 but want "chr" prepended to chromosome names, we can create
a function in the file `dm6.py` called `fasta_postprocess` that does this:

```python
def fasta_postprocess(origfn, newfn):
    shell("""zcat {origfn} | sed "s/>/>chr/g" | gzip -c > {newfn}  && rm {origfn}""")
```

Note that the file is uncompressed in order to fix the chromosome naming and
then recompressed into the final output. In this case the original file is
removed.

## Requirements

All files are required to be gzipped. Note that there is a `unzip_fasta` rule
that will temporarily unzip the fasta for things like bowtie2 index building
that don't support gzipped input.

Other than that, it's up to the user to decide what transformations (if any)
are required. Examples might include:

- exluding contigs
- removing or editing problematic genes that have transcripts on both strands
* renaming chromosomes (e.g., prepend "chr")
* remove unnecessary annotations (e.g., keep only cds/exon/transcript/gene features)
