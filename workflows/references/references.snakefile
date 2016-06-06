import os
import sys
import yaml
import importlib


def resolve_name(name):
    """
    Imports a specific object from a dotted path and returns just that object.

    From nose.utils.resolve_name (with the logging parts taken out) which in
    turn is from unittest.TestLoader.loadTestByName
    """
    parts = name.split('.')
    parts_copy = parts[:]
    while parts_copy:
        try:
            module = __import__('.'.join(parts_copy))
            break
        except ImportError:
            del parts_copy[-1]
            if not parts_copy:
                raise
    parts = parts[1:]
    obj = module
    for part in parts:
        obj = getattr(obj, part)
    return obj


def download_and_postprocess(assembly, field, outfile):
    """
    For the given assembly (dm6, hg19, etc) and field (gtf, fasta, etc) and
    desired output file, download the specified URL for the field to a temp
    file and apply the specified post-processing function (if any).

    The temp file is named after the output file; if the output file ends with
    .gz then the temp file will so as to not confuse other tools.

    This provides a mechanism for general downloading across assemblies while
    also allowing a sort of plugin customization option.
    """

    def default_postprocess(origfn, newfn):
        """
        If no other postprocess function is defined, then simply move the original
        to the new.
        """
        shell("mv {origfn} {newfn}")
    name = config[assembly][field].get('postprocess', None)
    if name is None:
        func = default_postprocess
    else:
        func = resolve_name(name)
    url = config[assembly][field]['url']
    if outfile.endswith('.gz'):
        gz = '.gz'
    else:
        gz = ''
    shell("wget {url} -O- > {outfile}.tmp{gz}")
    func(outfile + '.tmp' + gz, outfile)

ext_mapping = {
    'gtf': '.gtf.gz',
    'fasta': '.fa',
    'transcriptome': '.transcriptome.fa',
}

# Build the targets based on what's in the config file, rather than assume
# we want to build everything.
targets = []
for assembly, d in config.items():
    for key in d.keys():

        # TODO: eventually we'll want to add targets for DEXSeq, intergenic,
        # and astalavista variants
        if key == 'gtf':
            targets.append('{assembly}/{assembly}.gtf.gz'.format(assembly=assembly))

        if key == 'transcriptome':
            targets.append('{assembly}/{assembly}.transcriptome.fa.gz'.format(assembly=assembly))
            for index in d['transcriptome'].get('indexes', []):
                if index == 'kallisto':
                    targets.append('{assembly}/kallisto/{assembly}.transcriptome.idx'.format(assembly=assembly))
                else:
                    raise ValueError("transcriptome index %s not currently supported" % index)

        if key == 'fasta':
            targets.append(
                '{assembly}/{assembly}.fa.gz'.format(
                    assembly=assembly, ext=ext_mapping[key]))
            for index in d['fasta'].get('indexes', []):
                if index == 'bowtie2':
                    targets += expand('{assembly}/bowtie2/{assembly}.{n}.bt2', assembly=assembly, n=[1, 2, 3, 4])
                elif index == 'hisat2':
                    targets += expand('{assembly}/hisat2/{assembly}.{n}.ht2', assembly=assembly, n=range(1, 9))
                else:
                    raise ValueError("fasta index %s not currently supported" % index)


rule all:
    input: targets

# The downloading rules all have the same general form and support arbitrary
# post-processing functions to be specified in the config file.

rule download_fasta:
    output: '{assembly}/{assembly}.fa.gz'
    run:
        download_and_postprocess(wildcards.assembly, 'fasta', output[0])

rule unzip_fasta:
    input: rules.download_fasta.output
    output: temporary('{assembly}/{assembly}.fa')
    shell: 'gunzip -c {input} > {output}'

rule download_gtf:
    output: '{assembly}/{assembly}.gtf.gz'
    run:
        download_and_postprocess(wildcards.assembly, 'gtf', output[0])


rule download_transcriptome:
    output: '{assembly}/{assembly}.transcriptome.fa.gz'
    run:
        download_and_postprocess(wildcards.assembly, 'transcriptome', output[0])


rule bowtie_index:
    output: expand('{{assembly}}/bowtie2/{{assembly}}.{n}.bt2', n=[1, 2, 3, 4])
    input: rules.unzip_fasta.output
    log: '{assembly}/bowtie2/{assembly}.log'
    shell:
        '''
        bowtie2-build {input} {wildcards.assembly}/bowtie2/{wildcards.assembly} > {log} 2> {log}
        '''


rule hisat2_index:
    output: expand('{{assembly}}/hisat2/{{assembly}}.{n}.ht2', n=range(1, 9))
    log: '{assembly}/hisat2/{assembly}.log'
    input: rules.unzip_fasta.output
    shell:
        '''
        hisat2-build {input} {wildcards.assembly}/hisat2/{wildcards.assembly} > {log} 2> {log}
        '''


rule kallisto_index:
    output: '{assembly}/kallisto/{assembly}.transcriptome.idx'
    input: '{assembly}/{assembly}.transcriptome.fa.gz'
    log: '{assembly}/kallisto/{assembly}.log'
    shell:
        '''
        kallisto index -i {output} --make-unique {input} > {log} 2> {log}
        '''
# vim: ft=python
