#!/usr/bin/env bash

# Example bash script to install Pace to run bare-metal on ursa cluster.
# If you haven't already install miniconda for environment management (or use venv)

set -e -x

# module load necessary system software
module load python/3.11
module load openmpi
module load netcdf-c
module load cmake
module load cuda

# DO NOT load the gcc module, openmpi was compiled against the default gcc version,
# not the gcc in the module

# Assuming you want to use conda:
conda create -n pace python=3.11
conda activate pace

# # Else use a venv:
# python -m venv .vpace
# source .vpace/bin/activate

pip install --upgrade pip setuptools wheel

# If you want notebooks use these:
# pip install jupyter
# pip install ipyparallel

# Install cupy for GPUs
pip install cupy-cuda12x

# Install pace
cd pace
pip install .[test]
