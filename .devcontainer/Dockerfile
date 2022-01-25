FROM --platform=linux/x86_64 shinaoka/dcore_dev_triqs3

ENV PYTHONUNBUFFERED=1

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
RUN pip3 install scipy h5py!=2.10.0 toml dcorelib mpi4py matplotlib mypy pytest