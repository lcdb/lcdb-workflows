# Wrapper for Preseq c_curve

Preseq method that generates a complexity plot of the genome library based on observed data.

## Input
* Coordinate-sorted BAM. If paired-end, include `-P` in params.extra to run in
paired-end mode.

## Output
* two column text file displaying the *total reads* and a corresponding number of *distinct reads*

## Threads
Threads not supported.

## Example
```python
rule preseq_ccurve:
    input: '{sample}.sorted.bam'
    output: '{sample}.complexity_output.txt'
    wrapper:
        "file://path/to/preseq/observed_complexity"
```
