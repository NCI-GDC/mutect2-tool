.PHONY: docker-*
docker-login:
	@echo
	docker login -u="${QUAY_USERNAME}" -p="${QUAY_PASSWORD}" quay.io

.PHONY: build build-*
build: build-multi-mutect2

build-%:
	@echo
	@echo -- Building docker --
	@make -C $* build-docker NAME=$*

.PHONY: publish publish-% publish-release publish-release-%

publish-staging: publish-staging-multi-mutect2

publish-staging-%:
	@echo
	@make -C $* publish-staging

publish-release: publish-release-multi-mutect2

publish-release-%:
	@echo
	@make -C $* publish-release

