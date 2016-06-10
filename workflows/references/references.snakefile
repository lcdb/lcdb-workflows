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

d = {}
for i in config['references']:
    try:
        tag = i['tag']
    except KeyError:
        tag = 'default'
    d[(i['assembly'], tag)] = i


def download_and_postprocess(outfile):
    def default_postprocess(origfn, newfn):
        """
        If no other postprocess function is defined, then simply move the original
        to the new.
        """
        shell("mv {origfn} {newfn}")

    base = os.path.relpath(outfile, config['data_dir'])
    assembly = os.path.basename(os.path.dirname(base))
    basename = os.path.basename(outfile)

    tag = basename.split('_', 1)[-1].split('.')[0]
    block = d[(assembly, tag)]

    name = block.get('postprocess', None)
    if name is None:
        func = default_postprocess
    else:
        func = resolve_name(name)
    url = block['url']
    if outfile.endswith('.gz'):
        gz = '.gz'
    else:
        gz = ''
    shell("wget {url} -O- > {outfile}.tmp{gz}")
    func(outfile + '.tmp' + gz, outfile)


if not os.path.exists(config['data_dir']):
    os.makedirs(config['data_dir'])

workdir: config['data_dir']

# Build the targets based on what's in the config file, rather than assume
# we want to build everything.
ext_mapping = {
    'gtf': '.gtf.gz',
    'fasta': '.fa.gz',
}
index_mapping = {
    'bowtie2': ('{assembly}/bowtie2/{assembly}{tag}.{n}.bt2', dict(n=[1, 2, 3, 4])),
    'hisat2': ('{assembly}/hisat2/{assembly}{tag}.{n}.ht2', dict(n=range(1, 9))),
    'kallisto': ('{assembly}/kallisto/{assembly}{tag}.idx', dict()),
}
targets = []
for block in config['references']:
    try:
        tag = '_' + block['tag']
    except KeyError:
        tag = '_default'
    ext = ext_mapping[block['type']]
    assembly = block['assembly']
    targets.append('{assembly}/{assembly}{tag}{ext}'.format(**locals()))
    if block['type'] == 'fasta':
        indexes = block.get('indexes', [])
        for index in indexes:
            pattern, kwargs = index_mapping[index]
            kwargs = kwargs.copy()
            targets.extend(expand(pattern, assembly=assembly, tag=tag,  **kwargs))


rule all:
    input: targets

# The downloading rules all have the same general form and support arbitrary
# post-processing functions to be specified in the config file.

rule download_fasta:
    output: '{assembly}/{assembly}{tag}.fa.gz'
    run:
        download_and_postprocess(output[0])

rule unzip_fasta:
    input: rules.download_fasta.output
    output: temporary('{assembly}/{assembly}{tag}.fa')
    shell: 'gunzip -c {input} > {output}'

rule download_gtf:
    output: '{assembly}/{assembly}{tag}.gtf.gz'
    run:
        download_and_postprocess(output[0])


rule bowtie_index:
    output: expand('{{assembly}}/bowtie2/{{assembly}}{{tag}}.{n}.bt2', n=[1, 2, 3, 4])
    input: rules.unzip_fasta.output
    log: '{assembly}/bowtie2/{assembly}{tag}.log'
    shell:
        '''
        bowtie2-build {input} {wildcards.assembly}/bowtie2/{wildcards.assembly}{wildcards.tag} > {log} 2> {log}
        '''


rule hisat2_index:
    output: expand('{{assembly}}/hisat2/{{assembly}}{{tag}}.{n}.ht2', n=range(1, 9))
    log: '{assembly}/hisat2/{assembly}{tag}.log'
    input: rules.unzip_fasta.output
    shell:
        '''
        hisat2-build {input} {wildcards.assembly}/hisat2/{wildcards.assembly}{wildcards.tag} > {log} 2> {log}
        '''


rule kallisto_index:
    output: '{assembly}/kallisto/{assembly}{tag}.idx'
    input: '{assembly}/{assembly}{tag}.fa.gz'
    log: '{assembly}/kallisto/{assembly}{tag}.log'
    shell:
        '''
        kallisto index -i {output} --make-unique {input} > {log} 2> {log}
        '''
# vim: ft=python
