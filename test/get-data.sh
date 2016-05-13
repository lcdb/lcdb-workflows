#!/bin/bash
set -e

HERE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd $HERE
wget http://helix.nih.gov/~dalerr/lcdb-workflows-data/pasilla.tar
tar -xvf pasilla.tar
rm pasilla.tar
