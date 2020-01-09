FROM ubuntu:16.04
LABEL maintainer="sli6@uchicago.edu"
LABEL version="nightly-2016-02-25-gf39d340"
LABEL description="GATK3 nightly-2016-02-25-gf39d340"

RUN apt-get update -y && \
    apt-get install -y python-pip python-dev openjdk-8-jre-headless

ADD GenomeAnalysisTK.jar /opt/