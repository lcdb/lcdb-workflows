#!/bin/bash

# Downloads data used in the pasilla R package . . . but only the first 250k
# reads on chr2L.

set -e 
set -o pipefail
HERE=$(pwd)
CHROM="chr2L"
NREADS=250000
for i in \
    treated1 \
    treated2 \
    treated3 \
    untreated1 \
    untreated2 \
    untreated3 \
    untreated4 \
; do
    echo "${i}"
    mkdir -p "${i}"
    cd "${HERE}/${i}"
    if [ -e "${i}_R1.fastq" ]; then
        echo "${i}_R1.fastq exists; skipping"
    else
        samtools view -b "http://www-huber.embl.de/pub/DEXSeq/analysis/brooksetal/bam/${i}.bam" "$CHROM" \
            | bedtools sample -i stdin -n "$NREADS" \
            | bedtools bamtofastq -i stdin -fq "${i}.fastq.tmp" && mv "${i}.fastq.tmp" "${i}_R1.fastq" && rm ${i}.bam.bai
    fi
    cd ${HERE}
done

