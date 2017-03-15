FROM quay.io/jeremiahsavage/cdis_base

USER root
RUN apt-get update && apt-get install -y --force-yes \
    openjdk-8-jre-headless \
    libpq-dev \
    python-psycopg2

USER ubuntu
ENV HOME /home/ubuntu

ENV mutect2-tool 0.8e

RUN mkdir -p ${HOME}/tools/mutect2-tool
ADD docker/GenomeAnalysisTK.jar ${HOME}/tools/
ADD mutect2-tool ${HOME}/tools/mutect2-tool/
ADD setup.* ${HOME}/tools/mutect2-tool/
ADD requirements.txt ${HOME}/tools/mutect2-tool/

RUN /bin/bash -c "source ~/.virtualenvs/p3/bin/activate \
    && cd ~/tools/mutect2-tool \
    && pip install -r ./requirements.txt"

WORKDIR ${HOME}
