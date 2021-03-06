#!/bin/bash
set -e

HERE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd $HERE
if [ ! -e pasilla ]; then
    wget http://helix.nih.gov/~dalerr/lcdb-workflows-data/pasilla.tar
    tar -xvf pasilla.tar
    rm pasilla.tar
else
    echo "'pasilla' directory already exists"
fi
