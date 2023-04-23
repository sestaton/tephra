FROM ubuntu:latest

LABEL maintainer "S. Evan Staton"
LABEL image_type "Tephra: A tool for discovering transposable elements and describing patterns of genome evolution"

ENV LC_ALL="C"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    build-essential wget zlib1g-dev libgd-dev unzip libncurses5 libncurses5-dev libdb-dev git cpanminus \
    libexpat1 libexpat1-dev libxml2-dev libssl-dev openjdk-8-jre-headless gcc-multilib 

RUN echo "n" | cpanm -q -n App::cpm Data::Stag DB_File Bio::Root::Version Bio::SearchIO::blastxml Bio::SearchIO::hmmer 

RUN git clone https://github.com/sestaton/tephra.git \
    && cd tephra \
    && cpm install -g --show-build-log-on-failure \
    && perl Makefile.PL \
    && make install

RUN cpanm --uninstall HTTP::Cookies HTML::Tree HTTP::Negotiate WWW::RobotRules LWP \
    && apt-get remove -y git cpanminus \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf ~/.cpanm \	   
    && rm -rf ~/tephra
