FROM ubuntu:16.04

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
RUN apt-get update && apt-get install -y software-properties-common -y apt-transport-https
RUN apt-add-repository multiverse
RUN add-apt-repository restricted
RUN add-apt-repository -y ppa:webupd8team/java
RUN echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections
RUN apt-get update 
RUN apt-get install -y \
     wget \
     unzip \
     sudo \
     make \
     gcc \
     libdb-dev \
     libgd-dev \
     zlib1g-dev \
     libxml2-dev \
     libexpat1-dev \
     pkg-config \
     graphviz \
     libssl-dev \
     cpanminus




# Download and install

WORKDIR /opt/
COPY CRISPRCasMeta /opt/CRISPRCasMeta
WORKDIR /opt/CRISPRCasMeta/
RUN chmod 777 installer.sh
COPY CRISPRCasFinder.pl /opt/CRISPRCasMeta/
ENV PATH="/opt/CRISPRCasMeta:${PATH}"
RUN ./installer.sh

RUN echo "Complete 2"
WORKDIR /opt/CRISPRCasMeta/
