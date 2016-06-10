# Wrapper for Preseq gc_extrap

Preseq method that predicts the genomic coverage from deep sequencing based on observed data.

## Input
* BAM/BED: *Sorted* according to chromosome and start position

## Output
* future_coverage.txt: four column text file displaying the *total bases*, *expected covered bases* and the corresponding *lower/upper 95% confidence interval*

## Threads
Threads not supported.

## Params
* <code>-B, -bam</code>: should be included if input is a BAM file
* <code>-o, -output</code>: allows one to specify a name for the output file
* <code>-P, -pe</code>: if input is paired end preseq will register both mapped ends separately
* <code>-c, -cval</code>: confidence intervals (Default at 0.95)
* <code>-Q, -quick</code>: quick mode to estimate yield without bootstrapping for confidence intervals

## Example
<pre><code>
rule preseq_ccurve:
    input: 
		"mapped/{sample}.sorted.bam"
    output:
        "mapped/future_coverage.txt"
    params:
        "-o"
    wrapper:
        "file://path/to/preseq/estimated_coverage"
</code></pre>
