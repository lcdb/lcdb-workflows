# Wrapper for HISAT2



## Input
Index prefix

## Output
* List the output files in bullet form
* And talk about what they do, see exampe below:
* img: an image file showing the orientation of the workflow

## Threads
*Mention if multiple processors can be utilized.*  
*Provide the optimal number of processors to use.*

## Params
* Include the parameters in bullet form
* And describe why they will be needed, look at example below:
* <code>--force</code>: is always automatically specified but allows Snakemake to determine custom dependent files

## Example
Include the example Snakemake rule in <code>code</code> form, see example below:
<pre><code>
rule samtools_index:
    input: 
		"mapped/{sample}.sorted.bam"
    output:
        "mapped/{sample}.sorted.bam.bai"
    params:
        "optional params string"
    wrapper:
        "/path/to/wrapper/location"
</code></pre>
