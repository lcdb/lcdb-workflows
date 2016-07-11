#!/usr/bin/env python
"""
get data from geo
dm_rnaseq_se - GSE72947 - sir2
dm_rnaseq_pe - GSE49587 - smn
dm_chipseq_se - GSE49779 - tango
mm_rnaseq_se - GSE65976 - Maf1
hg_rnaseq_se - GSE73317 - brd4

"""
import os
import subprocess as sp
import string

'''
Requires fastq-dump. Loading with module as below didn't work.
sp.call(['module load sratoolkit'], shell=True)
adding sra-tools with conda '''

''' samples passed as dictionary'''

dm_rnaseq_se = {
        'WT1': 'SRR2352582',
        'WT2': 'SRR2352583',
        'KO1': 'SRR2352586',
        'KO2': 'SRR2352587',
    }

dm_rnaseq_pe = {
        'WT1': 'SRR948304',
        'WT2': 'SRR948305',
        'KO1': 'SRR948306',
        'KO2': 'SRR948307',
    }

dm_chipseq_se = {
        'untreated1': 'SRR1164481',
        'untreaded2': 'SRR1164482',
        'treated1': 'SRR1164484',
        'treated2': 'SRR1164485',
    }

mm_rnaseq_se = {
        'WT1': 'SRR1805875',
        'WT2': 'SRR1805876',
        'KO1': 'SRR1805878',
        'KO2': 'SRR1805879',
    }

hg_rnaseq_se = {
        'CON1': 'SRR2480700',
        'CON2': 'SRR2480701',
        'KD1': 'SRR2480704',
        'KD2': 'SRR2480705',
    }

libsize = {
        'small': '1000000',
        'large': '10000000',
        'full': 'full',
    }

def download(SPECIES, EXP, LAYOUT, samples, libsize):
    for size, READS in libsize.iteritems():
        f = open('sampletable.txt', 'wb')
        f.write("{}\t{}\t{}\n".format('sample_id', 'condition', 'replicate'))
    
        for sample_id, sample in samples.iteritems():
            print SPECIES, EXP, LAYOUT, sample, size
            if size == 'full':
                sp.call(['fastq-dump', '--gzip', '--split-files',  sample])
            else:
                sp.call(['fastq-dump', '--gzip', '-X', READS, '--split-files',  sample])
    
            dir = "-".join([SPECIES, EXP, LAYOUT, size])
            condition = sample_id.translate(None, string.digits)
            replicate = sample_id.translate(None, string.letters)
    
            try:
                os.stat(dir)
            except:
                os.makedirs(dir)
    
            if LAYOUT == 'SE':
                os.rename((sample + '_1.fastq.gz'),
                          os.path.join(dir, sample_id + '_R1.fastq.gz'))
            else:
                os.rename((sample + '_1.fastq.gz'),
                          os.path.join(dir, sample_id + '_R1.fastq.gz'))
                os.rename((sample + '_2.fastq.gz'),
                          os.path.join(dir, sample_id + '_R2.fastq.gz'))
            f.write("{}\t{}\t{}\n".format(sample_id, condition, replicate))
        f.close()
        os.rename('sampletable.txt', os.path.join(dir, 'sampletable.txt'))

download('dm', 'rnaseq', 'SE', dm_rnaseq_se, libsize)
download('dm', 'rnaseq', 'PE', dm_rnaseq_pe, libsize)
download('dm', 'chipseq', 'SE', dm_chipseq_se, libsize)
download('mm', 'rnaseq', 'SE', mm_rnaseq_se, libsize)
download('hg', 'rnaseq', 'SE', hg_rnaseq_se, libsize)
