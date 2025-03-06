ARG REGISTRY=docker.osdc.io/ncigdc
ARG BASE_CONTAINER_VERSION=latest

FROM quay.io/ncigdc/gatk:3.7 AS gatk
FROM ${REGISTRY}/python3.9-builder:${BASE_CONTAINER_VERSION} as builder
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

COPY ./ /mutect2_tool

WORKDIR /mutect2_tool

RUN pip install tox && tox -e build

FROM ${REGISTRY}/python3.9:${BASE_CONTAINER_VERSION}

LABEL org.opencontainers.image.title="mutect2_tool" \
      org.opencontainers.image.description="The GATK3 MuTect2 pipeline employs a 'Panel of Normals' to identify additional germline mutations. This panel is generated using TCGA blood normal genomes from thousands of individuals that were curated and confidently assessed to be cancer-free. This method allows for a higher level of confidence to be assigned to somatic variants that were called by the MuTect2 pipeline." \
      org.opencontainers.image.source="https://github.com/NCI-GDC/mutect2-tool" \
      org.opencontainers.image.vendor="NCI GDC"

COPY --from=builder /mutect2_tool/dist/*.whl /mutect2_tool/
COPY requirements.txt /mutect2_tool/

WORKDIR /mutect2_tool

RUN pip install --no-deps -r requirements.txt \
	&& pip install --no-deps *.whl \
	&& rm -f *.whl requirements.txt

USER app

ENTRYPOINT ["/tini", "--", "mutect2_tool"]

CMD ["--help"]
