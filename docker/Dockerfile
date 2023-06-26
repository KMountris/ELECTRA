FROM ubuntu:22.04

MAINTAINER Konstantinos A. Mountris <konstantinos.mountris@gmail.com>

# Copy src files from host to image.
COPY ./src /home/src

# Install basic libs
RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y apt-utils \
    git \
    g++ \
    make \
    wget \
    nano \
    xz-utils \
    mlocate \
    libssl-dev \
    libncurses5 \
    libncurses5-dev \
    libomp-dev \
    libgmp-dev \
    libmpfr-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libboost-thread-dev \
    libboost-program-options-dev \
    libopenblas-base \
    libopenblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    libsuitesparse-dev \
    libsuperlu-dev

# Download and install cmake
RUN cd /home \
    && wget https://github.com/Kitware/CMake/releases/download/v3.16.5/cmake-3.16.5.tar.gz \
    && tar xvf cmake-3.16.5.tar.gz \
    && rm cmake-3.16.5.tar.gz \
    && mv cmake-3.16.5 cmake \
    && mkdir /home/cmake/build \
    && cd /home/cmake/build \
    && ../bootstrap --parallel=4 \
    && make -j4 && make install \
    && cd ../.. \
    && rm -rf cmake

# Download and install tinyxml2
RUN cd /home \
    && git clone https://github.com/leethomason/tinyxml2.git \
    && mkdir /home/tinyxml2/build \
    && cd /home/tinyxml2/build \
    && cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=OFF\
    && make -j4 && make install \
    && cd ../.. \
    && rm -rf tinyxml2

# Download and install termcolor
RUN cd /home \
    && git clone https://github.com/ikalnytskyi/termcolor.git \
    && mkdir /home/termcolor/build \
    && cd /home/termcolor/build \
    && cmake .. -DCMAKE_BUILD_TYPE=RELEASE \
    && make -j4 && make install \
    && cd ../.. \
    && rm -rf termcolor

# Download and install Eigen
RUN cd /home \
    && wget https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz \
    && tar xvf eigen-3.3.9.tar.gz \
    && rm eigen-3.3.9.tar.gz \
    && mv eigen-3.3.9 eigen \
    && mkdir /home/eigen/build \
    && cd /home/eigen/build \
    && cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DEIGEN_TEST_CXX11=ON -DEIGEN_TEST_OPENMP=ON -DEIGEN_TEST_SSE2=ON -DEIGEN_TEST_SSE3=ON -DEIGEN_TEST_SSE4_1=ON -DEIGEN_TEST_SSE4_2=ON \
    && make -j4 && make install \
    && cd ../.. \
    && rm -rf eigen

# Download and install Armadillo
RUN cd /home \
    && wget http://sourceforge.net/projects/arma/files/armadillo-10.2.1.tar.xz \
    && tar xvf armadillo-10.2.1.tar.xz \
    && rm armadillo-10.2.1.tar.xz \
    && mv armadillo-10.2.1 armadillo \
    && mkdir /home/armadillo/build \
    && cd /home/armadillo/build \
    && cmake .. -DCMAKE_BUILD_TYPE=RELEASE \
    && make -j4 && make install \
    && cd ../.. \
    && rm -rf armadillo

# Download and install CGAL
RUN cd /home \
    && wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-5.0.1/CGAL-5.0.1.tar.xz \
    && tar xvf CGAL-5.0.1.tar.xz \
    && rm CGAL-5.0.1.tar.xz \
    && mv CGAL-5.0.1 cgal \
    && mkdir /home/cgal/build \
    && cd /home/cgal/build \
    && cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DCGAL_HEADER_ONLY=ON -DWITH_Eigen3=ON -DWITH_GMP=ON -DWITH_MPFR=ON \
    && make -j4 && make install \
    && cd ../.. \
    && rm -rf cgal

# Build install IMP
RUN cd /home \
    && mkdir /home/src/imp/build \
    && cd /home/src/imp/build \
    && cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_APPS=OFF -DBUILD_DOC=OFF -DBUILD_TESTS=OFF -DBUILD_SHARED_LIBS=OFF \
    && make -j4 && make install \
    && cd ../../..

# Build install CLOUDEA
RUN cd /home \
    && mkdir /home/src/cloudea/build \
    && cd /home/src/cloudea/build \
    && cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DCLOUDEA_USE_CGAL=ON -DBUILD_APPS=OFF -DBUILD_DOC=OFF -DBUILD_SHARED_LIBS=OFF \
    && make -j4 && make install \
    && cd ../../..

# Build install ELECTRA
RUN cd /home \
    && mkdir /home/electra \
    && cd /home/electra \
    && cmake ../src/electra -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_DOC=OFF -DELECTRA_WITH_CUDA=OFF -DBUILD_SHARED_LIBS=OFF \
    && make -j4 \
    && cd ../../..

# Delete src files
RUN rm -rf /home/src

# Add shortcut to ELECTRA applications
# RUN echo "alias ELECTRA-console='/home/electra/bin/ELECTRA-console'" >> ~/.bashrc
RUN echo "alias ElectraSim='/home/electra/bin/ElectraSim'" >> ~/.bashrc
RUN echo "alias ElectraPre='/home/electra/bin/ElectraPre'" >> ~/.bashrc

# Finished successfully
# constructing ELECTRA image.
CMD ["bash"]
