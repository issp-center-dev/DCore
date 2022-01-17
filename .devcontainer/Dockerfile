FROM python:3.9-slim

ENV PYTHONUNBUFFERED=1

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    curl \
    ca-certificates \
    git \
    zip \
    vim \
    gosu \
    libopenblas-base libopenblas-dev libhdf5-103 libhdf5-dev \
    cmake pkg-config gfortran \
    openmpi-bin libopenmpi-dev \
    && \
    apt-get clean && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* # clean up


RUN mkdir /tmp/dcore_work
WORKDIR /tmp/dcore_work
COPY ./ .

RUN pip3 install --upgrade pip setuptools
RUN pip3 install git+https://github.com/shinaoka/dcorelib.git
RUN python3 setup.py egg_info
RUN pip3 install `grep -v '^\[' src/dcore.egg-info/requires.txt`

# Create non-root user
ARG NB_USER=vscode
ARG NB_UID=1000
RUN useradd -u $NB_UID -m $NB_USER -s /bin/bash && \
    echo 'vscode ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER $NB_USER
ENV PATH "/home/${NB_USER}/.local/bin:${PATH}"
ENV PYTHONPATH "/home/${NB_USER}/work/src:${PYTONPATH}"

# for vscode
RUN mkdir /home/${NB_USER}/work

# For DCore
ENV DCORE_TRIQS_COMPAT 1
#test