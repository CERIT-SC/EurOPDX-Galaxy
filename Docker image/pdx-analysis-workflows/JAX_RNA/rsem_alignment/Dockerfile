FROM continuumio/miniconda:latest

WORKDIR /home/tools

COPY environment.yml ./
COPY rsem_alignment.py ./

RUN conda env create -f environment.yml
RUN rm environment.yml
RUN echo "source activate app-conda-env" > ~/.bashrc

ENV PATH /opt/conda/envs/app-conda-env/bin:$PATH