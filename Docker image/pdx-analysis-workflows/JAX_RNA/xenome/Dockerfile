FROM python:3

RUN apt -y update && apt-get -y install g++ cmake libboost-all-dev pandoc zlib1g-dev libbz2-dev libsqlite3-dev git

RUN git clone https://github.com/data61/gossamer.git
RUN cd gossamer
RUN mkdir build
WORKDIR /gossamer/build
RUN cmake ..
RUN make
RUN make install

RUN apt -y purge g++ cmake libboost-all-dev pandoc zlib1g-dev libbz2-dev libsqlite3-dev git
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /
COPY filter_rna_coverage.py ./
