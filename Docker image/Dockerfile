FROM ubuntu:20.04

RUN apt-get update
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
ENV TZ=Europe/London
RUN apt-get update
RUN apt-get install -y g++ cmake libboost-all-dev pandoc zlib1g-dev \
	libbz2-dev libsqlite3-dev apt-utils git python3 python3-distutils \
	python2 curl default-jre

# To install xenome from gossamer repo :
RUN git clone https://github.com/data61/gossamer.git
# this line causes trouble during installation on ubuntu > 16
RUN sed -i '98d' /gossamer/src/CMakeLists.txt
RUN cd gossamer
RUN mkdir build
WORKDIR /gossamer/build
RUN cmake -DBUILD_tests=OFF .. 
RUN make
RUN make install

WORKDIR /

RUN git clone --branch release_21.09 --depth 1 https://github.com/galaxyproject/galaxy

COPY ./galaxy-configs/install.sh galaxy/
RUN mkdir import/ 

WORKDIR galaxy

# pipeline tools and custom welcome.html
COPY ./galaxy-configs/tool_conf.xml ./config/tool_conf.xml
COPY ./galaxy-configs/galaxy.yml ./config/galaxy.yml
COPY ./galaxy-configs/dh-importer ./tools/dh-importer
COPY ./pdx-analysis-workflows ./tools/pdx-analysis-workflows
COPY ./galaxy-configs/static/ ./static/
COPY ./galaxy-configs/pipeline-architecture.svg ./static/images/

RUN chmod +x install.sh

RUN ./install.sh

EXPOSE 8080

CMD ["./run.sh"]

