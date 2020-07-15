FROM python:3.7-stretch
LABEL maintainer="sli6@uchicago.edu"
LABEL version="nightly-2016-02-25-gf39d340"
LABEL description="GATK3 nightly-2016-02-25-gf39d340"

RUN apt-get update -y \
    && apt-get install -y \
    openjdk-8-jre-headless \
    && apt-get clean

COPY . /opt/

WORKDIR /opt

RUN chmod +x multi_mutect2.py \
	&& ln -s /opt/multi_mutect2.py /bin/multi_mutect2

ENTRYPOINT ["/bin/multi_mutect2"]

CMD ["--help"]
