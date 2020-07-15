VERSION = 1.0
REPO = multi-mutect2

GIT_SHORT_HASH:=$(shell git rev-parse --short HEAD)
COMMIT_HASH:=$(shell git rev-parse HEAD)

DOCKER_URL := quay.io/ncigdc
DOCKER_IMAGE_COMMIT := ${DOCKER_URL}/${REPO}:${COMMIT_HASH}
DOCKER_IMAGE_LATEST := ${DOCKER_URL}/${REPO}:latest
DOCKER_IMAGE_STAGING := ${DOCKER_URL}/${REPO}:staging
DOCKER_IMAGE_PRODUCTION := ${DOCKER_URL}/${REPO}:${VERSION}

.PHONY: version version-* name
name:
	@echo ${NAME}

version:
	@echo --- VERSION: ${VERSION} ---

version-docker:
	@echo ${DOCKER_IMAGE_COMMIT}
	@echo ${DOCKER_IMAGE}

.PHONY: run
run:
	@docker run --rm ${DOCKER_IMAGE_LATEST}

.PHONY: build build-*

build: build-docker

build-docker:
	@echo
	@echo -- Building docker --
	docker build . \
		--file ./Dockerfile \
		--build-arg NAME=${NAME} \
		-t "${DOCKER_IMAGE_COMMIT}" \
		-t "${DOCKER_IMAGE_LATEST}"

.PHONY: docker-*
docker-login:
	@echo
	docker login -u=${QUAY_USERNAME} -p=${QUAY_PASSWORD} quay.io


.PHONY: publish-staging publish-release
publish: docker-login
	docker push ${DOCKER_IMAGE_COMMIT}

publish-staging: publish
	docker tag ${DOCKER_IMAGE_LATEST} ${DOCKER_IMAGE_STAGING}
	docker push ${DOCKER_IMAGE_STAGING}

publish-release: publish
	docker tag ${DOCKER_IMAGE_LATEST} ${DOCKER_IMAGE_PRODUCTION}
	docker push ${DOCKER_IMAGE_LATEST}
	docker push ${DOCKER_IMAGE_PRODUCTION}

