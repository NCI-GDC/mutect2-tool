FROM quay.io/jeremiahsavage/cdis_base

USER root
RUN apt-get update && apt-get install -y --force-yes \
    openjdk-8-jre-headless

USER ubuntu
ENV HOME /home/ubuntu

ENV mutect-tool 0.4c

RUN mkdir -p ${HOME}/tools/mutect-tool
ADD docker/GenomeAnalysisTK.jar ${HOME}/tools/
ADD mutect-tool ${HOME}/tools/mutect-tool/
ADD setup.* ${HOME}/tools/mutect-tool/

RUN /bin/bash -c "source ${HOME}/.local/bin/virtualenvwrapper.sh \
    && source ~/.virtualenvs/p3/bin/activate \
    && cd ~/tools/mutect-tool \
    && pip install -e ."

WORKDIR ${HOME}
