FROM flatironinstitute/triqs

# Copy host directory
RUN echo $HOME
COPY ./ /home/triqs/src/DCore
RUN sudo chown -R triqs /home/triqs/src/DCore
RUN ls -la /home/triqs/src/DCore

# Install poetry for build
ENV PATH $PATH:/home/triqs/.local/bin
RUN curl -kL https://bootstrap.pypa.io/get-pip.py | python3
RUN pip3 install sphinx wild_sphinx_theme matplotlib --user

#RUN mkdir /home/triqs/dcore_build
WORKDIR /home/triqs/src/DCore
RUN ls -la
RUN python3 setup.py install --user
RUN sphinx-build -b html doc html
RUN find /home/triqs/src/DCore/html

#RUN mkdir $HOME/build/DCore
#RUN mkdir $HOME/build/DCore
#WORKDIR $HOME/build/DCore
#ENV CTEST_OUTPUT_ON_FAILURE=1
#RUN cmake -DCMAKE_VERBOSE_MAKEFILE=ON $HOME/src/DCore -DTRIQS_PATH=$HOME/opt/triqs -DCMAKE_INSTALL_PREFIX=$HOME/opt/DCore -DBUILD_DOC=ON \
 #&& make -j 2 all && make -j 2 install
#RUN find .
