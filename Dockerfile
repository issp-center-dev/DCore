FROM shinaoka/triqs1.4_dmft_alps

# Copy host directory
COPY ./ $HOME/src/DCore

WORKDIR $HOME/src/DCore
RUN mkdir $HOME/build/DCore
WORKDIR $HOME/build/DCore
ENV CTEST_OUTPUT_ON_FAILURE=1
RUN cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DTRIQS_PATH=$HOME/opt/triqs $HOME/src/DCore -DCMAKE_INSTALL_PREFIX=$HOME/opt/DCore -DBUILD_DOC=ON\
 && make -j 2 install
