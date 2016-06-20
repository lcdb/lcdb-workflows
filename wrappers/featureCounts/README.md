# Wrapper for featureCounts

Read summarizaton/quantification is required for a number of downstream 
analyses such as gene expression and histone modification.  The output 
is a counts table where the number of reads assigned to each feature in
each library is recorded. *FeatureCounts* can be used to perform this task
much faster than HTSeq without compromising on accuracy. 

[Link to homepage](http://bioinf.wehi.edu.au/featureCounts/)

[Link to manual](http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf)

## Input
* `BAM`: BAM file; paired-end reads will be automatically re-ordered location sorted
* `GTF`: Genomic Feature annotations are used to identify gene locations

## Output
* `Counts`: *counts.txt* contains the number of reads assigned to each meta-feature

## Threads
Multiple processors supported.

## Params
* `extra`: a string that will be passes verbatim to featureCounts.

## Example

```python
rule featureCounts:
    input:
		annoGTF="{sample}.gtf",
		mappingFile="{sample}.bam"
    output:
        counts="{sample}_counts.txt"
	threads: 5
    params: 
		extra="-t exon -g gene_id "
    wrapper:
		"/path/to/wrapper/location"
```
