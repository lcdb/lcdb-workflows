#!/usr/bin/env python

import argparse
import subprocess as sp
import os
import uuid
import sys

HERE = os.path.realpath(os.path.dirname(__file__))
ENV_NAME = 'lcdb-workflows-%s-env' % os.environ['USER']
REQUIREMENTS = os.path.join(HERE, 'requirements.txt')
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'

ap = argparse.ArgumentParser()

ap.add_argument(
    'workflow',
    default='test/Snakefile',
    help='Which workflow to test (path to snakefile)')

ap.add_argument(
    '--build-env', action='store_true',
    help='Build environment "%s"' % ENV_NAME)

ap.add_argument(
    '--no-test', action='store_true',
    help="Make the temp bash script, but don't actually run the test. "
    "Useful in combination with --build-env")

ap.add_argument(
    '--clean', action='store_true',
    help="Run the 'clean' rule of the snakefile before calling it a second time")

ap.add_argument(
    '--sbatch', action='store_true',
    help='Print the sbatch command required to run the tests')

ap.add_argument(
    '--threads',
    type=int,
    help='''Threads to run. If --sbatch is specified, then request this many
    --cpus-per-task; if $SLURM_CPUS_PER_TASK exists (that is, you're on an
    interactive node) use that; otherwise assume running on a node with this
    many cores.  Snakemake will be run with the -j argument set to this many
    threads. If --cluster is specified, this has no effect (the cluster
    submission command limits the number of jobs to 999)''')

ap.add_argument(
    '--mem',
    default='32g',
    help='Memory to request. Only has an effect if --sbatch also specified.')

ap.add_argument(
    '--cluster',
    action='store_true',
    help='Instead of submitting a single job to the cluster, build a submit '
    'script and cluster config and submit that to the cluster')

ap.add_argument(
    '--config',
    help='''Specify a default config file. The default is to use the
    config.yaml file in the same directory as the specified workflow.''')

args = ap.parse_args()

if args.config:
    CONFIG = args.config

else:
    CONFIG = os.path.join(os.path.dirname(args.workflow), 'config.yaml')

# If no threads supplied and we're on an interactive node, use as many CPUs as
# are allocated.
if args.threads is None:
    try:
        threads = os.environ['SLURM_CPUS_PER_TASK']
    except KeyError:
        threads = 1
else:
    threads = args.threads


def env_exists(name):
    """
    See if the environment already exists by checking output of `conda env
    list`.
    """
    return any(
        map(
            lambda x: name in x,
            sp.check_output(
                ['conda', 'env', 'list'], universal_newlines=True
            ).splitlines()
        )
    )

def create_env(name):
    """
    Remove environment if it exists and build a new one.
    """
    if env_exists(name):
        print(
            YELLOW
            + 'conda env '
            + name + ' exists, removing and rebuilding based on '
            + requirements(repo)
            +  ENDC)
        sp.check_call([
            'conda', 'remove', '-n', name, '--all'])
    sp.check_call([
        'conda', 'create', '-n', name, '--file', requirements(repo), '-c', 'bioconda', '-c', 'r', 'python=3', '-y'])


if args.build_env:
    create_env(ENV_NAME)

DIRECTORY = os.path.dirname(args.workflow)
SNAKEMAKE = "snakemake --directory {DIRECTORY} -s {args.workflow}".format(**locals())

if args.clean:
    CLEAN = SNAKEMAKE + ' clean'
else:
    CLEAN = ""

# This will be saved as a temp file and either run or printed out in the sbatch
# command.
script = """\
#!/bin/bash
set -eo pipefail
source activate {ENV_NAME}
{CLEAN}
{SNAKEMAKE} -j {threads} -p -r -T --verbose --configfile {CONFIG}
"""
script_name = 'lcdb-workflows-submit-%s.sh' % (str(uuid.uuid4()).split('-')[0])
with open(script_name, 'w') as fout:
    fout.write(script.format(**locals()))
print(BLUE + "Wrote " + script_name  + ENDC)

if args.sbatch:
    cmd = [
        'sbatch',
        '--mem', args.mem,
        '--cpus-per-task', str(threads),
        '--mail-type=END,FAIL',
        script_name]
    sp.check_call(cmd)
    sys.exit(0)

elif args.cluster:
    # run the separate submit script so that it can do the config parsing and
    # cluster-config.yaml building work
    cmd = [
        'python', 'lcdb/lcdb-submit.py',
        '--config', CONFIG,
        '--workflow', args.workflow,
        '--cluster-config', 'test-cluster-config.yaml',
        '--output', script_name,
        '--pre-block', CLEAN
    ]
    sp.check_call(cmd)

if not args.no_test:
    print(BLUE + "Running " + script_name + ENDC)
    sp.check_call(['bash', script_name])
