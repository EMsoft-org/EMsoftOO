FROM ubuntu:jammy

ARG TARGETARCH 
ARG DEBIAN_FRONTEND=noninteractive

RUN if [ "$TARGETARCH" = "arm64" ]; then \
 apt-get update && apt-get -y upgrade \
 && apt-get install -y \
 git \
 vim \
 wget \
 cmake \
 gcc \
 gfortran \
 ninja-build \
 build-essential \
 libopenblas-dev \
 apt-utils \
 opencl-headers \
 ocl-icd-opencl-dev \
 ocl-icd-libopencl1 \
 clinfo ; \
 elif [ "$TARGETARCH" = "amd64" ]; then \
 apt-get update && apt-get -y upgrade \
 && apt-get install -y \
 git \
 vim \
 wget \
 cmake \
 gcc \
 gfortran \
 ninja-build \
 build-essential \
 libopenblas-dev \
 apt-utils \
 opencl-headers \
 ocl-icd-opencl-dev \
 ocl-icd-libopencl1 \
 libpocl-dev \ 
 intel-opencl-icd \
 clinfo \
; fi

# intel-opencl-icd makes iGPUs work from Intel (tested on i5-8400)
# libpocl-dev is for CPU on both Intel and AMD

# set some environmental variables for the Nvidia container
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility
# set up OpenCL for Nvidia
RUN mkdir -p /etc/OpenCL/vendors && \
    echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# clone EMsoft and set up SDK Debug/Release
RUN mkdir /home/EMs \
 && cd /home/EMs \
 && git clone --branch developOO https://github.com/EMsoft-org/EMsoftSuperbuild.git \
 && git clone https://github.com/EMsoft-org/EMsoftData.git \
 && git clone https://github.com/EMsoft-org/EMsoftOO.git \
 && mv EMsoftSuperbuild EMsoftOOSuperbuild \
 && mkdir EMsoftOOBuild && mkdir GenEMsoftData && cd EMsoftOOSuperbuild && mkdir Debug Release

# EMsoftSuperbuild
RUN cd /home/EMs/EMsoftOOSuperbuild/Debug/ \
 && cmake -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DCMAKE_BUILD_TYPE=Debug ../ -G Ninja && ninja \
 && cd ../Release \
 && cmake -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DCMAKE_BUILD_TYPE=Release ../ -G Ninja && ninja

# EMsoftOOBuild ... must use dynamic "shared" libraries for f90wrap later on

RUN if [ "$TARGETARCH" = "arm64" ]; then \
 cd /home/EMs/EMsoftOOBuild/ && mkdir Debug Release && cd Debug \
 && cmake -DCMAKE_BUILD_TYPE=Debug -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DBUILD_SHARED_LIBS=ON \
 ../../EMsoftOO -G Ninja \
 && ninja \
 && cd ../Release \
 && cmake -DCMAKE_BUILD_TYPE=Release -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DBUILD_SHARED_LIBS=ON \
  ../../EMsoftOO -G Ninja \
 && ninja \
 && cd /home/EMs/GenEMsoftData && mkdir EMXtal \
 && cd ../../ && mkdir .config && cd .config && mkdir EMsoft && cd EMsoft && mkdir tmp && cd ../../ ; \
 elif [ "$TARGETARCH" = "amd64" ]; then \
 cd /home/EMs/EMsoftOOBuild/ && mkdir Debug Release && cd Debug \
 && cmake -DCMAKE_BUILD_TYPE=Debug -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DBUILD_SHARED_LIBS=ON \
 ../../EMsoftOO -G Ninja \
 && ninja \
 && cd ../Release \
 && cmake -DCMAKE_BUILD_TYPE=Release -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DBUILD_SHARED_LIBS=ON \
  ../../EMsoftOO -G Ninja \
 && ninja \
 && cd /home/EMs/GenEMsoftData && mkdir EMXtal \
 && cd ../../ && mkdir .config && cd .config && mkdir EMsoft && cd EMsoft && mkdir tmp && cd ../../ \
 ; fi

# add release version to path
ENV PATH /home/EMs/EMsoftOOBuild/Release/Bin:$PATH

# run terminal for user at /home/
WORKDIR /home/
CMD ["/bin/bash"]
