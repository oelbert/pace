#!/usr/bin/env bash

# Example bash script to install Pace to run bare-metal on Gaea's c4 cluster

set -e -x

# module load necessary system software
module load gcc/8
module load openmpi/gcc/4.1.6
module load boost/1.76.0
module load anaconda3/2024.10

# export CC=`which gcc`
# export CXX=`which g++`
# export MPICC=`which mpicc`
# export MPICXX=`which mpicxx`
# export DACE_compiler_cpu_executable=${CXX}
# export GT4PY_EXTRA_COMPILE_ARGS="-O3 -ffast-math -fprefetch-loop-arrays -funroll-all-loops"
# export OPENMP_CPPFLAGS="-fopenmp -fopenmp-simd"
# export OPENMP_LDFLAGS="-fopenmp -fopenmp-simd"

# create a conda environment for pace
conda create -y --name pace python=3.11.7

# enter the environment and update it
conda activate my_name
pip install --upgrade pip setuptools wheel

# install the Pace dependencies, GT4Py, and Pace
pip install -r requirements_dev.txt

# If you want to run notebooks:
pip install ipyparallel
