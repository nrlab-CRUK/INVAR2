
version ?= latest

build:
	docker build --tag crukcibioinformatics/invar2:${version} .

release: build
	docker push crukcibioinformatics/invar2:${version}

