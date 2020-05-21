FROM ubuntu:18.04
LABEL maintainer="sli6@uchicago.edu"
LABEL version="nightly-2016-02-25-gf39d340"
LABEL description="GATK3 nightly-2016-02-25-gf39d340"

RUN apt update -y \
    && apt install -y \
    openjdk-8-jre-headless \
    software-properties-common

RUN add-apt-repository ppa:deadsnakes/ppa \
    && apt update -y \
    && apt install python3.7 -y

ADD GenomeAnalysisTK.jar /opt/
ADD multi_mutect2_p3.py /opt/