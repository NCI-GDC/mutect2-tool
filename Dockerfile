FROM quay.io/ncigdc/gatk:3.7 AS gatk
FROM quay.io/ncigdc/python37 AS python
FROM openjdk:8u282-slim

ENV BINARY mutect2-tool
LABEL maintainer="sli6@uchicago.edu"
LABEL version="nightly-2016-02-25-gf39d340"
LABEL description="GATK3 nightly-2016-02-25-gf39d340"

COPY --from=python / /
COPY --from=gatk /usr/local/bin/ /usr/local/bin/

COPY ./dist/ /opt
WORKDIR /opt

RUN apt-get update \
	&& apt-get install make \
	&& rm -rf /var/lib/apt/lists/*

RUN make init-pip \
  && ln -s /opt/bin/${BINARY} /usr/local/bin/${BINARY}

ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
ENTRYPOINT ["/tini", "--", "mutect2_tool"]
CMD ["--help"]
