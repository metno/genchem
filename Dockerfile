FROM ubuntu:18.04

##########################################################################################
# Root is the default user of the docker container, which may pose a security risk.      #
# Please ensure this is acceptable before proceeding.                                    #
##########################################################################################

# Install and configure OS
RUN apt-get update
RUN apt-get install -y build-essential \
    apt-utils \
    wget \
    libgl1-mesa-glx \
    gfortran \
    linux-headers-generic \
    cmake \
    vim \
    unzip \
    pkgconf \
    libpng-dev \
    libfreetype6-dev \
    libfontconfig1 \
    libxrender1 \
    xauth \
    git \
    subversion

### Set up directory structure
RUN mkdir -p /Code/
RUN mkdir -p /Code/downloads/

### Download anaconda for Python
WORKDIR /Code/downloads/
RUN wget -q https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh

### Install anaconda
RUN bash Anaconda3-2018.12-Linux-x86_64.sh -b -p /Code/anaconda3
RUN rm Anaconda3-2018.12-Linux-x86_64.sh

### Configure system for anaconda
ENV PATH="/Code/anaconda3/bin:${PATH}"
RUN conda update -y conda

### Clone genchem from github (creates genchem directory)
WORKDIR /Code/
Run git clone https://github.com/metno/genchem

### Test box setup from https://genchem.readthedocs.io/en/latest/GenChemDoc_quickstart.html
WORKDIR /Code/genchem/box/
RUN ./scripts/box_setup.sh tmp_work
WORKDIR /Code/genchem/box/tmp_work
RUN ./do.testChems EmChem19a
