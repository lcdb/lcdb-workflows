FROM ubuntu:14.04

MAINTAINER Ryan Dale <dalerr@niddk.nih.gov>

#RUN apt-get update && apt-get install -y \
    build-essential \
    bzip2 \
    ca-certificates \
    git \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxtst6 \
    libxrender1 \
    mysql-client \
    wget \
    x11-apps \
    zlib1g-dev

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --no-check-certificate --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH /opt/conda/bin:$PATH
ENV DISPLAY :99
ENV LANG C.UTF-8
RUN mkdir /opt/lcdb
COPY requirements.txt /tmp/
RUN conda install -y -c bioconda -c r --file /tmp/requirements.txt
RUN apt-get install -y libxi6 xvfb xauth
RUN apt-get install -y libxext-dev libxrender-dev libxtst-dev
ADD xvfb_init /etc/init.d/xvfb
RUN chmod a+x /etc/init.d/xvfb
CMD ["/etc/init.d/xvfb", "start]
WORKDIR /opt/lcdb
