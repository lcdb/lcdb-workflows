# Wrapper for *addtoolnamehere*

*Include a brief description about the tool here.*

Input
=====
*Include what the input files will be and how the tool handles them.*

Output
======
* List the output files in bullet form
* And talk about what they do, see exampe below:
* img: an image file showing the orientation of the workflow

Threads
=======
*Mention if multiple processors can be utilized.*  
*Provide the optimal number of processors to use.*

Params
======
* Include the parameters in bullet form
* And describe why they will be needed, look at example below:
* <code>--force</code>: is always automatically specified but allows Snakemake to determine custom dependent files

Example
=======
Include the example Snakemake rule in <code>code</code> form, see example below:
<pre><code>
rule multiqc
    input: expand('samples/{sample}{suffix}', sample=samples, suffix=['_fastqc.zip', '.cutadapt.log'])
    output:
        html='multiqc_report.html'
    params:
        analysis_directory='samples'
        extra="--config=multiqc_config.yaml"
    wrapper:
        "file://path/to/multiqc"
</code></pre>
