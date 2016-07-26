import os
import pandas
import yaml
from jsonschema import validate, ValidationError
from snakemake.shell import shell
from snakemake.io import expand

ALIGNER_TAGS = {
    'bowtie2': 'bt2',
    'hisat2': 'ht2',
}

def validate_config(config, schema):
    schema = yaml.load(open(schema))
    cfg = yaml.load(open(config))
    try:
        validate(cfg, schema)
    except ValidationError as e:
        msg = '\nPlease fix %s: %s\n' % (config, e.message)
        raise ValidationError(msg)


def build_wrapper_for(source_dir, wrappers_dir):
    """
    Returns a `wrapper_for` function to be used in a workflow.

    Parameters
    ----------
    :source_dir: str
        Directory of the calling snakemake workflow. Typically this is obtained
        with the srcdir() built-in.
    :wrappers_dir: str
        Directory of wrappers relative to source dir
    """
    def wrapper_for(tool):
        return os.path.join(source_dir, wrappers_dir, tool)
    return wrapper_for


def build_params_for(config):
    """
    Returns a `params_for` function to be used in a workflow.

    Parameters
    ----------
    :config: dict
        The global config dictionary from a workflow
    """
    def params_for(rule, key):
        return  config.get('rules', {}).get(rule, {}).get('params', {}).get(key, '')
    return params_for


def build_threads_for(config):
    """
    Returns a `threads_for` function to be used in a workflow.

    Parameters
    ----------
    :config: dict
        The global config dictionary from a workflow
    """
    def threads_for(rule):
        return  config.get('rules', {}).get(rule, {}).get('threads', None)
    return threads_for


def workflow_helper_functions(config, source_dir, wrappers_dir):
    """
    One-stop-shop for building helper functions.

    Parameters
    ----------
    :config: dict
        The global config dictionary from a workflow
    :source_dir: str
        Directory of the calling snakemake workflow. Typically this is obtained
        with the srcdir() built-in.
    :wrappers_dir: str
        Directory of wrappers relative to source dir

    Returns
    -------
    wrappers_for, params_for, and threads_for functions.
    """
    return (
        build_wrapper_for(source_dir, wrappers_dir),
        build_params_for(config),
        build_threads_for(config),
    )


def load_sampletable(filename):
    """
    Load sampletable.

    TODO: validation will go here.
    """
    return pandas.read_table(filename, index_col=0)


def rscript(string, scriptname, log=None):
    """
    Saves the string as `scriptname` and then runs it

    Parameters
    ----------
    string : str
        Filled-in template to be written as R script

    scriptname : str
        File to save script to

    log : str
        File to redirect stdout and stderr to. If None, no redirection occurs.
    """
    with open(scriptname, 'w') as fout:
        fout.write(string)
    if log:
        _log = '> {0} 2>&1'.format(log)
    else:
        _log = ""
    shell('Rscript {scriptname} {_log}')


def aligner_tag(config):
    return ALIGNER_TAGS[config['rules']['align']['aligner']]


def aligner(config):
    return config['rules']['align']['aligner']


def aligner_index(config):
    aligner = config['rules']['align']['aligner']
    if aligner == 'hisat2':
        return expand('{data_dir}/{prefix}.{n}.{tag}',
                      prefix=config['rules']['align']['prefix'],
                      n=range(1, 9),
                      tag=aligner_tag(config),
                      data_dir=config['data_dir'])

    elif aligner == 'bowtie2':
        return expand('{data_dir}/{prefix}.{n}.{tag}',
                      prefix=config['rules']['align']['prefix'],
                      n=range(1, 5),
                      tag=aligner_tag(config),
                      data_dir=config['data_dir'])

    else:
        raise ValueError('Unsupported aligner "{}"'.format(aligner))


def fastq_screen_config_inputs(config):
     """
     Returns indexes and a table (as string) that can be written directly to
     a config file.

     Use indexes as the input to a rule; use the table to build a config file.
     """
     indexes = config['rules']['fastq_screen']['indexes']
     aligner = config['rules']['fastq_screen']['aligner']
     inputs = []
     for label, prefix in indexes.items():
         prefix = os.path.join(config['data_dir'], prefix)
         inputs.extend(expand(prefix + '.{n}.bt2', n=[1, 2, 3, 4]))
     return inputs
