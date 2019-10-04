FROM ubuntu:18.04

LABEL maintainer "S. Evan Staton"
LABEL image_type "Tephra: A tool for discovering transposable elements and describing patterns of genome evolution"

ARG LC_ALL="C"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    build-essential wget default-jre zlib1g-dev libgd-dev unzip libncurses5 libncurses5-dev libdb-dev git \
    cpanminus libexpat1 libexpat1-dev libidn11 perl-doc \
    && echo "n" | cpanm -q -n Data::Stag DB_File Bio::Root::Version Bio::SearchIO::blastxml Bio::SearchIO::hmmer \
    && git clone https://github.com/sestaton/tephra.git \
    && cd tephra \
    && cpanm -q --installdeps . \
    && perl Makefile.PL \
    && make install \
    && cpanm --uninstall HTTP::Daemon HTTP::Cookies HTML::Tree HTTP::Negotiate WWW::RobotRules LWP \
    && apt-get remove -y git cpanminus unzip \
    && rm -rf /var/lib/apt/lists/* \
    && cd .. && rm -rf tephra
