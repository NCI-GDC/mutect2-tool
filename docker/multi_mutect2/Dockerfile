FROM ubuntu:16.04
LABEL maintainer="sli6@uchicago.edu"
LABEL version="1.4"
LABEL description="Python wrapper of multithreading `MuTect2` in GATK3 nightly-2016-02-25-gf39d340."

RUN apt-get update -y \
    && apt-get install -y \
    python-dev \
    openjdk-8-jre-headless

ADD GenomeAnalysisTK.jar /opt/
ADD multi_mutect2.py /opt/