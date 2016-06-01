# Key Features for RSeQC
Comprehensively evaluate high throughput sequence data, especially RNA-seq data

1. *bam_stat.py* checks the mapping statistics of reads that:
    * QC fail
	* Uniquely map
	* Splice map
	* Map in a proper pair
2. *inner_distance.py* estimates the inner distance distribution between paired reads
    * Estimated inner distance should be consistent with gel size selection
	* Detects structural variation or aberrant splicing in RNA-seq data
3. *geneBody_coverage.py* scales all transcripts to 100 nt and calculates the number of reads covering each nucleotide position
    * Generates a plot illustrating the coverage profile along the gene body  
![](https://github.com/lcdb/lcdb-workflows/wrappers/rseqc/plotA.png)
4. *read_distribution.py* 

## Basic RSeQC Modules
* Sequence quality
* Nucleotide composition bias
* PCR bias
* GC bias

## RNA-seq Specific
* Sequencing saturation
* Mapped reads distribution
* Coverage uniformity
* Strand specificity
* Transcript level RNA itegrity

## FastQC vs RSeQC
* FastQC only focuses on raw sequence-related metrics
    + Not enough to ensure the usability of RNA-seq data
* RSeQC accesses the quality of RNA-seq experiments
