FROM marcdegraef/emsoftoo_sdk:buildx-latest

ARG TARGETARCH 
ARG DEBIAN_FRONTEND=noninteractive

# clone EMsoft and set up SDK Debug/Release
RUN cd ~/EMs \
 && git clone https://github.com/EMsoft-org/EMsoftData.git \
 && git clone https://github.com/EMsoft-org/EMsoftOO.git \
 && mkdir EMPlay && mkdir EMsoftOOBuild && mkdir EMXtal

RUN cd ~/EMs/EMsoftOOBuild/ && mkdir Debug Release && cd Debug \
 && cmake -DCMAKE_BUILD_TYPE=Debug -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DBUILD_SHARED_LIBS=OFF \
 ../../EMsoft -G Ninja \
 && ninja \
 && cd ../Release \
 && cmake -DCMAKE_BUILD_TYPE=Release -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DBUILD_SHARED_LIBS=OFF \
  ../../EMsoft -G Ninja \
 && ninja
 
# add release version to path
ENV PATH ~/EMs/EMsoftOOBuild/Release/Bin:$PATH
# add backup path of EMsoft
ENV EMSOFTPATHNAME ~/EMs/EMsoftOO

# run terminal for user at /home/
WORKDIR /home/
CMD ["/bin/bash"]
