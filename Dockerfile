FROM python:3.12-slim-bookworm@sha256:28cf028e5a544e92dbe11450debd93dd5eb70eaf3179a9e878cfaee426556b3b

RUN apt-get update && apt-get install -y make \
    software-properties-common \
    libgeos-dev \
    libopenmpi3 \
    libopenmpi-dev \
    libhdf5-serial-dev \
    netcdf-bin \
    libnetcdf-dev \
    libproj-dev \
    proj-bin \
    proj-data \
    python3 \
    python3-pip \
    git

RUN python -m pip install --upgrade pip setuptools wheel

ENV SETUPTOOLS_SCM_PRETEND_VERSION_FOR_NDSL=2025.10.00

COPY . /pace

RUN cd /pace && \
    python -m pip install -e .[dev]

RUN cd / && \
    git clone https://github.com/ai2cm/fv3net

ENV CFLAGS="-I/usr/include -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H=1"


RUN python -m ensurepip --upgrade && \
    python -m pip install \
    matplotlib==3.10.0 \
    ipyparallel==8.4.1 \
    jupyterlab==3.4.4 \
    cartopy==0.23.0 \
    jupyterlab_code_formatter==1.5.2 \
    isort==5.10.1 \
    black==22.3.0 \
    /fv3net/external/vcm

RUN python -m pip install pybind11==2.13.6

ENV PYTHONPATH=/fv3net/external/fv3viz:/pace/external/gt4py/src

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
