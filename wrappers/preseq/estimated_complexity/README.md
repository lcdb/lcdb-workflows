# Wrapper for Preseq lc_extrap

Preseq method used to estimate the future yield of the genome library based on observed data.

## Input
* BAM/BED: *Sorted* according to chromosome and start position

## Output
* future_yield.txt: four column text file displaying the *total reads*, the *expected distinct* reads and its corresponding *lower/upper 95% Confidence Interval*

## Threads
Threads not supported.

## Params
* <code>-B, -bam</code>: should be included if input is a BAM file
* <code>-o, -output</code>: allows one to specify a name for the output file
* <code>-P, -pe</code>: if input is paired end preseq will register both mapped ends separately

## Example
<pre><code>
rule preseq_lcextrap:
    input: 
		"mapped/{sample}.sorted.bam"
    output:
        "mapped/future_yield.txt"
    params:
        "-B -o future_yield.txt"
    wrapper:
        "file://path/to/preseq/estimated_complexity"
</code></pre>