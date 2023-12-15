FROM marcdegraef/emsoftoo_sdk:buildx-latest

ARG TARGETARCH 
ARG DEBIAN_FRONTEND=noninteractive

# clone EMsoft and set up SDK Debug/Release
RUN cd /home/EMs \
 && git clone https://github.com/EMsoft-org/EMsoftData.git \
 && git clone https://github.com/ZacharyVarley/EMsoftOO.git \
 &&  mkdir EMsoftOOBuild 

RUN cd /home/EMs/EMsoftOOBuild/ && mkdir Debug Release && cd Debug \
 && cmake -DCMAKE_BUILD_TYPE=Debug -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DBUILD_SHARED_LIBS=OFF \
 ../../EMsoftOO -G Ninja \
 && ninja \
 && cd ../Release \
 && cmake -DCMAKE_BUILD_TYPE=Release -DEMsoftOO_SDK=/opt/EMsoftOO_SDK -DBUILD_SHARED_LIBS=OFF \
  ../../EMsoftOO -G Ninja \
 && ninja
 
# add release version to path
ENV PATH /home/EMs/EMsoftOOBuild/Release/Bin:$PATH
# add backup path of EMsoft
ENV EMSOFTPATHNAME /home/EMs/EMsoftOO

# install a new user
ARG user=EMuser
ARG uid=501
ARG gid=3000

# Add new user with our credentials
ENV USERNAME ${user}
RUN useradd -m $USERNAME && \
        echo "$USERNAME:$USERNAME" | chpasswd && \
        usermod --shell /bin/bash $USERNAME && \
        usermod  --uid ${uid} $USERNAME && \
        groupmod --gid ${gid} $USERNAME

USER ${user}

# create EMuser workfolders and set up EMsoftConfig.json file 
RUN mkdir /home/${user}/XtalFolder \
 && mkdir /home/${user}/EMPlay \
 && mkdir /home/${user}/.config \
 && mkdir /home/${user}/.config/EMsoft \
 && mkdir /home/${user}/.config/EMsoft/tmp \
 && cp /home/EMs/EMsoftOO/EMsoftDockerConfig.template /home/${user}/.config/EMsoft/EMsoftConfig.json

# run terminal for EMuser at /home/${user}
WORKDIR /home/${user}
CMD ["/bin/bash"]
