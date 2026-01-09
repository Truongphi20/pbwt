FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    autotools-dev \
    libtool \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    bzip2 \
    make \
    sudo \
    git \
    wget \
    ca-certificates \
    libboost-program-options-dev \
    libboost-iostreams-dev \
    libboost-serialization-dev \
    && rm -rf /var/lib/apt/lists/*

# Build HTSlib
WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -xjf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    ./configure --enable-libcurl && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd .. && rm -rf htslib-1.16*

# Install macs
RUN git clone https://github.com/gchen98/macs.git && \
    cd macs && \
    make && \
    chmod 777 macs && \
    cp macs /usr/bin && \
    cd .. && rm -rf macs

# Install cmake
RUN mkdir -p /usr/local
RUN cd /usr/local &&\
    wget https://github.com/Kitware/CMake/releases/download/v4.0.2/cmake-4.0.2-linux-x86_64.sh &&\
    chmod +x cmake-4.0.2-linux-x86_64.sh &&\
    ./cmake-4.0.2-linux-x86_64.sh --skip-license


# Install required packages to debug
RUN apt-get update && apt-get install -y --no-install-recommends openssh-client \
                                                gdb \
                                                bcftools \
                                                less \
                                                bash-completion

# Install Python
RUN sudo apt-get update && sudo apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    python-is-python3 \
    && sudo rm -rf /var/lib/apt/lists/*


# Create vscode user
RUN useradd -m -s /bin/bash vscode && \
    echo "vscode ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/vscode && \
    chmod 0440 /etc/sudoers.d/vscode

USER vscode
WORKDIR /workspace
