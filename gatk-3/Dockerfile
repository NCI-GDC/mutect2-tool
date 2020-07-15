FROM ubuntu:18.04
LABEL maintainer="sli6@uchicago.edu"
LABEL version="nightly-2016-02-25-gf39d340"
LABEL description="GATK3 nightly-2016-02-25-gf39d340"

RUN apt-get update -y \
    && apt-get install -y \
    openjdk-8-jre-headless

COPY GenomeAnalysisTK.jar /opt

ENTRYPOINT ["java"]

CMD ["-jar", "/opt/GenomeAnalysisTK.jar"]
