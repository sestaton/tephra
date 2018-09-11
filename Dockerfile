FROM ubuntu:18.04

LABEL maintainer "S. Evan Staton"
LABEL image_type "Tephra: A tool for discovering transposable elements and describing patterns of genome evolution"

RUN apt-get update \
    && apt-get upgrade -y -qq \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    build-essential wget default-jre zlib1g-dev unzip libncurses5 libncurses5-dev libdb-dev git cpanminus libexpat1 libexpat1-dev libidn11 \
    && cpanm Data::Stag DB_File \
    && echo "n" | cpanm -q -n Bio::Root::Version \
    && git clone https://github.com/sestaton/tephra.git \
    && cd tephra \
    && cpanm -q --installdeps . \
    && perl Makefile.PL \
    && make install \
    && wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz \
    && tar xzf muscle3.8.31_i86linux64.tar.gz \
    && mv muscle3.8.31_i86linux64 /usr/local/bin/muscle \
    && cpanm --uninstall HTTP::Daemon HTTP::Cookies HTML::Tree HTTP::Negotiate WWW::RobotRules LWP \
    && apt-get remove -y git cpanminus unzip \
    && rm -rf /var/lib/apt/lists/* \
    && cd .. && rm -rf tephra
