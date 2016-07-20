#!/usr/bin/env python
'''
Read a sampletable and download data from sra
By default downloads the full set
If libsize is specified then download full and a sebset specified by libsize
Need fastq_dump which can be installed using bionconda sra-tools
'''

import sys
import os
import string
import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description='Download fastq files from SRA.')
parser.add_argument('--libsize', default='full', help='library size (default: full)')
args = parser.parse_args()

fns = os.listdir("sampletable-dir")
#dir = [i.replace(".txt","") for i in fns]

if args.libsize != 'full':
    lines = str(int(args.libsize)*4)

for i in fns:
    with open(os.path.join('sampletable-dir',i)) as f:
        dir = i.replace(".txt","")
        next(f)
        rows = ( line.split('\t') for line in f )
        d = { row[0]:row[1:] for row in rows }
    for sra_id, sample_id in d.iteritems():
        print sra_id, sample_id[0]
        sp.call(['fastq-dump', '--gzip', '--split-files',  sra_id])

        dir_full = ''.join([dir, '-full'])
        dir_subset = ''.join([dir, '-', args.libsize])

        if not os.path.exists(dir_full):
            os.makedirs(dir_full)
            
        if "-se" in dir:
            os.rename((sra_id + '_1.fastq.gz'),
            os.path.join(dir_full, sample_id[0] + '_R1.fastq.gz'))
            if args.libsize != 'full':
                if not os.path.exists(dir_subset):
                    os.makedirs(dir_subset)
                sp.check_call(' '.join(['gzip', '-cd',
                    os.path.join(dir_full, sample_id[0] +  '_R1.fastq.gz'),
                    '|', 'head', '-n', lines, '|', 'gzip', '>',
                    dir_subset + '/' + sample_id[0] +  '_R1.fastq.gz']), shell=True)
        else:
            os.rename((sra_id + '_1.fastq.gz'),
              os.path.join(dir_full, sample_id[0] + '_R1.fastq.gz'))
            os.rename((sra_id + '_2.fastq.gz'),
              os.path.join(dir_full, sample_id[0] + '_R2.fastq.gz'))
            if args.libsize != 'full':
                if not os.path.exists(dir_subset):
                    os.makedirs(dir_subset)
                sp.check_call(' '.join(['gzip', '-cd',
                    os.path.join(dir_full, sample_id[0] +  '_R1.fastq.gz'),
                    '|', 'head', '-n', lines, '|', 'gzip', '>',
                    dir_subset + '/' + sample_id[0] +  '_R1.fastq.gz']), shell=True)
                sp.check_call(' '.join(['gzip', '-cd',
                    os.path.join(dir_full, sample_id[0] +  '_R2.fastq.gz'),
                    '|', 'head', '-n', lines, '|', 'gzip', '>',
                    dir_subset + '/' + sample_id[0] +  '_R2.fastq.gz']), shell=True)
