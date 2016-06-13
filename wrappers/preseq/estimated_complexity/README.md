# Wrapper for Preseq lc_extrap

Preseq method used to estimate the future yield of the genome library based on observed data.

## Input
* Coordinate-sorted BAM. If paired-end, include `-P` in params.extra to run in
paired-end mode.
## Output

* future_yield.txt: four column text file displaying the *total reads*, the
  *expected distinct* reads and its corresponding *lower/upper 95% Confidence
  Interval*

## Threads
Threads not supported.

## Params
extra: a string that will be passed verbatim to gc_extrap

## Example
```python
rule preseq_lcextrap:
    input: '{sample}.sorted.bam'
    output: '{sample}/future_yield.txt'
    wrapper:
        "file://path/to/preseq/estimated_complexity"
```
