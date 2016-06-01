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
![](https://cloud.githubusercontent.com/assets/11708268/15725008/bce81f50-2817-11e6-9b03-f205bde446f3.png)
4. *read_distribution.py* calculates the fraction of reads mapped to the following regions:
    * Coding regions
	* 5' UTR exons
	* 3' URT exons
	* Introns
	* Intergenic regions *based on the gene model provided*
	    + Background noise levels can be estimated by adding custom gene models
5. *RPKM_saturation.py* determines the precision of estimated *Reads Per Kilobase of transcript per Million* (RPKM)
    * Estimation is performed at current sequencing depth by resampling the total mapped reads
	* Percent relative error **(100 * |RPKM~obs~ - RPKM~real~ | / RPKM~real~)** is used to measure RPKM  
![](https://cloud.githubusercontent.com/assets/11708268/15725749/4ba206ae-281b-11e6-9065-8264178f0aad.png)

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
