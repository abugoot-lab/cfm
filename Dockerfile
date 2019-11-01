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
RUN ./installer_UBUNTU.sh

RUN wget https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=CRISPRCasFinder.zip -O CRISPRCasFinder.zip
RUN unzip CRISPRCasFinder.zip
WORKDIR /opt/CRISPRCasFinder/
RUN chmod 777 installer_UBUNTU.sh
RUN ./installer_UBUNTU.sh

#COPY CRISPRCasFinder_modified.pl /opt/CRISPRCasFinder/
WORKDIR /opt/CRISPRCasMeta/
#ENV PATH="/opt/CRISPRCasFinder/bin:${PATH}"
