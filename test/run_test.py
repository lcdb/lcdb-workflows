#!/usr/bin/env python

import argparse
import subprocess as sp
import os
import uuid
import sys

env_name = 'lcdb-workflows-%s-env' % os.environ['USER']


class bcolors:
    "Fancy terminal color output. Thanks http://stackoverflow.com/a/287944"
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


ap = argparse.ArgumentParser()

ap.add_argument('repo', help='Path to lcdb-workflows directory to test')

ap.add_argument(
    '--build-env', action='store_true',
    help='Build environment "%s"' % env_name)

ap.add_argument(
    '--no-test', action='store_true',
    help="Make the temp bash script, but don't actually run the test. "
    "Useful in combination with --build-env")

ap.add_argument(
    '--clean', action='store_true',
    help="Delete the test data directory and re-download data.")

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
    threads.''')

ap.add_argument(
    '--mem',
    default='32g',
    help='Memory to request. Only has an effect if --sbatch also specified.')

ap.add_argument(
    '--workflow',
    default='test',
    help='Which workflow to test')

args = ap.parse_args()

# If no threads supplied and we're on an interactive node, use as many CPUs as
# are allocated.
if args.threads is None:
    try:
        threads = os.environ['SLURM_CPUS_PER_TASK']
    except KeyError:
        threads = 1
else:
    threads = args.threads


REPO = os.path.abspath(args.repo)

print(bcolors.OKBLUE + "Using repo: " + REPO + bcolors.ENDC)

def requirements(repo):
    "Path to requirements.txt in the repo"
    return os.path.join(repo, 'test', 'requirements.txt')

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

def create_env(name, repo):
    """
    Remove environment if it exists and build a new one.
    """
    if env_exists(name):
        print(
            bcolors.WARNING
            + 'conda env '
            + name + ' exists, removing and rebuilding based on '
            + requirements(repo)
            +  bcolors.ENDC)
        sp.check_call([
            'conda', 'remove', '-n', name, '--all'])
    sp.check_call([
        'conda', 'create', '-n', name, '--file', requirements(repo), '-c', 'bioconda', '-c', 'r', 'python=3', '-y'])


if args.build_env:
    create_env(env_name, REPO)

if args.clean:
    CLEAN = "snakemake clean --configfile config.yaml"
else:
    CLEAN = ""

# This will be saved as a temp file and either run or printed out in the sbatch
# command.
script = """\
#!/bin/bash

set -e

source activate {env_name}
cd {repo}/{args.workflow}
./get-data.sh
{CLEAN}
snakemake --unlock --configfile config.yaml
snakemake -j {threads} -pr --configfile config.yaml
"""
script_name = 'lcdb-workflows-submit-%s.sh' % (str(uuid.uuid4()).split('-')[0])
with open(script_name, 'w') as fout:
    fout.write(script.format(args=args, env_name=env_name, repo=REPO, threads=threads, CLEAN=CLEAN))
print(bcolors.OKBLUE + "Wrote " + script_name  + bcolors.ENDC)

if args.sbatch:
    cmd = [
        'sbatch',
        '--mem', args.mem,
        '--cpus-per-task', str(threads),
        '--mail-type=END,FAIL',
        script_name]
    sp.check_call(cmd)
    sys.exit(0)

elif not args.no_test:
    print(bcolors.OKBLUE + "Running " + script_name + bcolors.ENDC)
    sp.check_call(['bash', script_name])
