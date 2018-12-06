FROM shinaoka/dev_ubuntu_user

# TRIQS

# This forces an update of the cache when the master branch moves
ADD https://api.github.com/repos/TRIQS/triqs/git/refs/heads/1.4.x version.json

# Clone and build triqs/master
RUN git clone -b 1.4.x --single-branch https://github.com/TRIQS/triqs.git $HOME/src/triqs
RUN mkdir $HOME/build && mkdir $HOME/build/triqs
WORKDIR $HOME/build/triqs
RUN pwd
RUN cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/triqs -DCMAKE_VERBOSE_MAKEFILE=ON -DBuild_Documentation=ON -DDocWithCpp2doc=OFF $HOME/src/triqs
RUN make -j 2 install && make test && make clean
