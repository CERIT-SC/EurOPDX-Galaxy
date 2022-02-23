FROM continuumio/miniconda:latest

WORKDIR /

COPY environment.yml ./
COPY dh-importer.xml ./
COPY dh-importer.py ./
COPY create_tsv_output.py ./
COPY create_tsv_output.xml ./
COPY join_tsv_output.py ./
COPY join_tsv_output.xml ./

RUN conda env create -f environment.yml
RUN rm environment.yml
RUN echo "source activate app-conda-env" > ~/.bashrc

ENV PATH /opt/conda/envs/app-conda-env/bin:$PATH
