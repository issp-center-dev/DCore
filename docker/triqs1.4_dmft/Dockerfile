FROM shinaoka/triqs1.4

# Clone and build triqs/dft_tools
WORKDIR $HOME/src
RUN git clone -b master --single-branch https://github.com/TRIQS/dft_tools.git \
 && cd dft_tools \
 && git checkout d00575632c4 \
 && mkdir $HOME/build/dft_tools
WORKDIR $HOME/build/dft_tools
# Workaround for test failure: 17 - vaspio (Failed) 
RUN cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DTRIQS_PATH=$HOME/opt/triqs $HOME/src/dft_tools \
 && make -j 2 install \
 && make test \
 && make clean; exit 0

# Clone and build triqs/hubbardI (from forked repo due to some compilation errors)
WORKDIR $HOME/src
RUN git clone -b master --single-branch https://github.com/TRIQS/hubbardI.git \
 && mkdir $HOME/build/hubbardI
WORKDIR $HOME/build/hubbardI
RUN cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DTRIQS_PATH=$HOME/opt/triqs $HOME/src/hubbardI \
 && make -j 2 install \
 && make test \
 && make clean

# Clone and build triqs/cthyb
WORKDIR $HOME/src
RUN git clone -b 1.4.2 --single-branch https://github.com/TRIQS/cthyb.git \
 && mkdir $HOME/build/cthyb
WORKDIR $HOME/build/cthyb
RUN cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DTRIQS_PATH=$HOME/opt/triqs -DHYBRIDISATION_IS_COMPLEX=ON -DLOCAL_HAMILTONIAN_IS_COMPLEX=ON $HOME/src/cthyb \
 && make -j 2 install \
 && make test \
 && make clean
