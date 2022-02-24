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
    hdf5-tools \
    libopenblas-base \
    libopenblas-dev \
    libhdf5-103 \
    libhdf5-dev \
    libeigen3-dev \
    cmake \
    pkg-config \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \
    sudo \
    libboost-dev \
    libboost-mpi-dev \
    && \
    apt-get clean && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* # clean up


# pomerol
WORKDIR /root
RUN git clone https://github.com/aeantipov/pomerol.git && cd pomerol && git checkout 7a45b6a8b8254dcbf8 && cd .. && \
    mkdir pomerol.build && cd pomerol.build && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/pomerol -DPOMEROL_COMPLEX_MATRIX_ELEMENTS=ON ../pomerol && \
    make install && cd .. && rm -rf pomerol.build pomerol

# pomerol2dcore
WORKDIR /root
RUN git clone https://github.com/j-otsuki/pomerol2dcore.git && \
    mkdir pomerol2dcore.build && cd pomerol2dcore.build && \
    cmake -Dpomerol_DIR=/opt/pomerol -DCMAKE_INSTALL_PREFIX=/opt/pomerol2dcore ../pomerol2dcore && \
    make install

# Create non-root user
ARG NB_USER=vscode
ARG NB_UID=1000
RUN useradd -u $NB_UID -m $NB_USER -s /bin/bash && \
    echo 'vscode ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER $NB_USER
ENV PATH "/home/${NB_USER}/.local/bin:/opt/pomerol2dcore/bin:${PATH}"
ENV PYTHONPATH "/home/${NB_USER}/work/src:${PYTONPATH}"
ENV OMPI_MCA_btl_vader_single_copy_mechanism "none"

# for vscode
RUN mkdir /home/${NB_USER}/work

# For DCore
ENV DCORE_TRIQS_COMPAT 1

RUN pip3 install -U pip && \
    pip3 install scipy h5py toml tomli dcorelib mpi4py matplotlib pytest mypy sympy
