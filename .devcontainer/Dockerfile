FROM --platform=linux/x86_64 shinaoka/dcore_dev_triqs3

ENV PYTHONUNBUFFERED=1

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    gnuplot \
    ghostscript \
    && \
    apt-get clean && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* # clean up

# Create non-root user
ARG NB_USER=vscode
ARG NB_UID=1000
RUN useradd -u $NB_UID -m $NB_USER -s /bin/bash && \
    echo 'vscode ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER $NB_USER
ENV PATH "/home/${NB_USER}/.local/bin:/opt/pomerol2dcore/bin:${PATH}"
ENV PYTHONPATH "/home/${NB_USER}/work/src:${PYTONPATH}"

# for vscode
RUN mkdir /home/${NB_USER}/work

# For DCore
ENV DCORE_TRIQS_COMPAT 1

RUN pip3 install --upgrade pip
RUN pip3 install numpy scipy h5py!=2.10.0 toml dcorelib mpi4py matplotlib mypy pytest -U
RUN pip3 install sphinx wild_sphinx_theme versioneer -U

RUN echo "source /opt/triqs/share/triqsvars.sh" >> /home/${NB_USER}/.bashrc