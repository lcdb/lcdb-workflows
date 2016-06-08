# Preseq

*preseq* is an empirical Bayes method for predicting the 
molecular complexity of sequencing libraries or samples
on the basis of data from very shallow sequencing runs.

*In english*: provide preseq with either a **BED** or **BAM**
format and it will give you the following information:
<ul>
<li>yield estimates from experiment and subsampled experiments</li>
<li>predicted yield estimates from additional sequencing</li>
<li>predict the genomic coverage given current reads</li>
</ul>

## Basic Usage

*sort*: preseq requires a sorted BED or BAM file (recommend BED)
<pre><code>$ sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 input.bed > input.sort.bed</code></pre>

*c_curve*: generates a complexity plot of the genomic library
in BED or BAM format
<pre><code>$ ./preseq c_curve -o complexity_output.txt input.bed</code></pre>

*lc_extrap*: estimate the future yield of a genomic library using an initial experiment
<pre><code>$ ./preseq lc_extrap -o future_yield.txt input.bed</code></pre>

*bam2mr*: converts bam (sorted) to nr (recommended before running *gc_extrap*)
<pre><code>$ ./bam2mr sorted.bam > sorted.nr</code></pre>

*gc_extrap*: predicts the genomic coverage from deep sequencing based on the initial sample
<pre><code>$ ./preseq gc_extrap -o future_coverage.txt sorted.nr</code></pre>


**Note**: for additional information please refer to the [manual](http://smithlabresearch.org/wp-content/uploads/manual.pdf)

## Example

<pre><code>
rule preseq:
input: 
bed='samples/{sample}.bed'
output:
c_curve='samples/{sample}_complexity_output.txt'
lc_curve='samples/{sample}_future_yield.txt'
gc_extrap='samples/{sample}_future_coverage.txt'

</code></pre>

Branch Testing
