FROM --platform=linux/x86_64 shinaoka/dcore_dev_triqs3

ENV PYTHONUNBUFFERED=1


# Install requirement
RUN mkdir /tmp/setup
WORKDIR /tmp/setup
COPY ./ .

RUN python3 setup.py egg_info
RUN pip3 install `grep -v '^\[' src/*.egg-info/requires.txt`

# Create non-root user
ARG NB_USER=vscode
ARG NB_UID=1000
RUN useradd -u $NB_UID -m $NB_USER -s /bin/bash && \
    echo 'vscode ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER $NB_USER
ENV PATH "/home/${NB_USER}/.local/bin:/opt/pomerol2dcore/bin:${PATH}"
#ENV PYTHONPATH "/home/${NB_USER}/work/src:${PYTONPATH}"

# for vscode
RUN mkdir /home/${NB_USER}/work

# For DCore
#ENV DCORE_TRIQS_COMPAT 0

RUN echo "source /opt/triqs/share/triqsvars.sh" >> /home/${NB_USER}/.bashrc
