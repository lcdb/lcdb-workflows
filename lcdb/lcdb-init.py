#!/usr/bin/env python

import os
import logging
import argparse
import subprocess as sp
import sys

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(os.path.basename(__file__))

HERE = os.path.dirname(os.path.realpath(__file__))
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'


def run(cmds, **kwargs):
    """
    Print commands in blue before running them. kwargs are passed to
    subprocess.check_call.
    """
    logger.info(BLUE + ' '.join(cmds) + ENDC)
    sp.check_call(cmds, **kwargs)


def clone(dest, label, repo=None):
    if repo is None:
        repo = 'git@github.com:lcdb/lcdb-workflows.git'
    if os.path.exists(dest):
        print(RED + dest + ' already exists, aborting' + ENDC)
        sys.exit(1)
    run(['git', 'clone', repo, dest])
    run(['git', 'checkout', '-b', label], cwd=dest)
    generate_config(dest)
    run(['git', 'add', 'config.yaml'], cwd=dest)
    run(['git', 'commit', '-m', '"initialization of %s"' % label], cwd=dest)


def build_env(repo, label):
    run([

        'conda', 'create', '-n', label,
        '--file', os.path.join(repo, 'test', 'requirements.txt'),
        '-c', 'bioconda', '-c', 'r', 'python=3'
    ])


def generate_config(dest):
    """
    Sets up configuration. Eventually this will build a config file based on
    the schema.
    """
    # TODO: lots more to add here; currently just creates something to commit.
    run(['touch', 'config.yaml'], cwd=dest)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '--dest', default='.',
        help='Destination directory in which to clone workflows')
    ap.add_argument(
        '--label',
        help='Label for experiment. Will be used for git branch and conda '
        'environment.')
    ap.add_argument(
        '--no-build-env', action='store_true',
        help="Don't build a new environment")
    args = ap.parse_args()
    clone(args.dest, args.label)
    if not args.no_build_env:
        build_env(args.dest, args.label)
