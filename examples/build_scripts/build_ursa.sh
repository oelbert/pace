#!/usr/bin/env bash

# Example bash script to install Pace to run bare-metal on ursa cluster

set -e -x

# module load necessary system software
module load gcc
module load python/3.11
module load openmpi
module load netcdf-c
module load cmake
module load cuda


python -m venv .vpace
source .vpace/bin/activate

pip install --upgrade pip setuptools wheel

# If you want notebooks use these:
# pip install jupyter
# pip install ipyparallel

# Install cupy for GPUs
pip install cupy-cuda12x

# Install pace
cd pace
pip install .[test]

#If this fails, be sure that /usr/include/ is in your PATH