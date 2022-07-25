#!/bin/bash

docker run --rm -v $PWD:/home/jovyan/work -p 8888:8888 jupyter/scipy-notebook
