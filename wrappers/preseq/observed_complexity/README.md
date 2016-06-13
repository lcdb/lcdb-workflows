# Wrapper for Preseq c_curve

Preseq method that generates a complexity plot of the genome library based on observed data.

## Input
* BAM/BED: *Sorted* according to chromosome and start position

## Output
* complexity_output.txt: two column text file displaying the *total reads* and a corresponding number of *distinct reads*

## Threads
Threads not supported.

## Params
* <code>-B, -bam</code>: should be included if input is a BAM file
* <code>-o, -output</code>: allows one to specify a name for the output file
* <code>-P, -pe</code>: if input is paired end preseq will register both mapped ends separately

## Example
<pre><code>
rule preseq_ccurve:
    input: 
		sortBAM='{sample}.sorted.bam'
    output:
        complexOut='complexity_output.txt'
    wrapper:
        "file://path/to/preseq/observed_complexity"
</code></pre>
