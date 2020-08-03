FROM continuumio/miniconda:4.7.12

WORKDIR /repliseq-nf

COPY environment.yml /repliseq-nf/environment.yml

RUN apt-get update \
    && apt-get install -y procps \
    && apt-get clean -y \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda env create --name repliseq-nf -f environment.yml \
    && rm -rf /opt/conda/pkgs/*

ENV PATH /opt/conda/envs/repliseq-nf/bin:$PATH
