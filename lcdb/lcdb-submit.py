#!/usr/bin/env python
import os
import argparse
import yaml

ap = argparse.ArgumentParser()
ap.add_argument('--config', help='Global config file')
ap.add_argument('--cluster-config', help='Specify the generated cluster '
                'config file. Default is %(default)s.', default='cluster-config.yaml')
ap.add_argument('--workflow', help='Workflow to submit (path to snakefile)')
ap.add_argument('--output', help='Output file for generated script. Default is stdout')
ap.add_argument('--pre-block', 
                help="""String for the 'PRE_BLOCK' placeholder in the generated
                script, which is right before the snakemake call. Primarily
                used by the testing framework, for example to include
                a "snakemake clean" command prior to running.""",
                default="")
ap.add_argument('remainder', nargs=argparse.REMAINDER,
                help='Additional args are passed to snakemake')
args = ap.parse_args()

config = yaml.load(open(args.config))

# extract the cluster info.

cluster_config = {'__default__': config['cluster_default']}
for rule, block in config.get('rules', {}).items():
    if 'cluster' in block:
        cluster_config[rule] = block['cluster']

yaml.dump(cluster_config, open(args.cluster_config, 'w'), default_flow_style=False)

dirname = os.path.dirname(args.workflow)
REMAINDER = " ".join(args.remainder)
SNAKEFILE = args.workflow

wrapper_template = """\
#!/bin/bash

set -eo pipefail

source activate lcdb-workflows-$USER-env
{args.pre_block}
time snakemake --cluster-config {args.cluster_config} \\
    --cluster "sbatch {{cluster.args}} --cpus-per-task={{threads}}" \\
    --jobname "s.{{rulename}}.{{jobid}}.sh" \\
    -j 999 \\
    --rerun-incomplete \\
    -T \\
    --verbose \\
    --directory {dirname} \\
    -s {SNAKEFILE} \\
    {REMAINDER} \\
    > {args.workflow}.log 2>&1
""".format(**locals())

if not args.output:
    print(wrapper_template)
else:
    with open(args.output, 'w') as fout:
        fout.write(wrapper_template)
