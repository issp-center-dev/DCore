FROM flatironinstitute/triqs

RUN sudo apt update && sudo apt install vim make --yes

# Copy host directory
RUN echo $HOME
COPY ./ /home/triqs/src/DCore
RUN sudo chown -R triqs /home/triqs/src
RUN ls -la /home/triqs/src

# Install TRIQS/HubbardI
WORKDIR /home/triqs/src
RUN git clone https://github.com/TRIQS/hubbardI.git
RUN mkdir hubbardI.build
WORKDIR /home/triqs/src/hubbardI.build
RUN cmake -DCMAKE_CXX_COMPILER=mpicxx ../hubbardI
RUN make && sudo make install

# Install pip & some python libraries
ENV PATH $PATH:/home/triqs/.local/bin
RUN curl -kL https://bootstrap.pypa.io/get-pip.py | python3
RUN pip3 install sphinx wild_sphinx_theme matplotlib pytest numpy -U --user

# Install DCore
WORKDIR /home/triqs/src/DCore
RUN ls -la
RUN python3 setup.py install --user

# Build doc
RUN sphinx-build -b html doc html
RUN find /home/triqs/src/DCore/html

# Run tests
RUN pytest tests/non-mpi/*/*.py
RUN pytest tests/mpi/*/*.py
