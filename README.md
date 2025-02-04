# EurOPDX-Galaxy
This repository contains recipe for docker version of [Galaxy platform](www.galaxyproject.org), expanded with tools and workflows for PDX(Patient Derived Xenograft) model molecular data analysis.

To use this docker image, you need to have [Docker](https://docs.docker.com/get-docker/) and [Docker Compose](https://docs.docker.com/compose/install/) installed.

# Quick start

This is a tutorial for UNIX-like environment (Lunix, MacOS)

### Setup a work-dir
```bash
mkdir ~/docker-galaxy
cd ~/docker-galaxy
mkdir reference-data
```

### Download the docker-compose.yml

In the `~/docker-galaxy` directory 

```bash
wget https://raw.githubusercontent.com/BorisYourich/EurOPDX-Galaxy/main/docker-compose.yml
```

### Configure



```yaml
version: '3.5'

services:
  galaxy:
    container_name: galaxy
    restart: always
    image: edirex/galaxy:latest
    ports:
      - 8080:8080

    volumes:
      - /your/import/dir:/import
      - /your/dir/to/store/reference/files:/galaxy/reference-data/
```


If you do not wish to build the image yourself, which is unnecessary unless you need to change something radical in the image, then you only need to download the [docker-compose.yml](https://github.com/BorisYourich/EurOPDX-Galaxy/blob/main/docker-compose.yml) file and set the correct paths in the "volumes:". 
First line defines which directory will be used as the library import directory by the container i.e. path to the folder where the data you want to proces are located. Second line tells the container where it should store the reference data e.g. Genome indexes, created by the mapper tools, so that you do not need to create them multiple times, because everytime you stop the running container, it deletes it's stored data. Do not change the second part of the argument after the ":" as it is associated with the internal structure of the image.

After these modifications you are ready to start the container. In the terminal, move to the directory where you have downloaded the docker-compose.yml file and run "docker-compose up" command. You should see the image being downloaded, this will take some time, but it is only necessary to do it once. After download, the container will start and you can see the EurOPDX Galaxy instance on "localhost:8080" on your machine.
