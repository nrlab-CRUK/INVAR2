
version ?= latest

build:
	docker build --tag nrlabcruk/invar2:${version} .

release: build
	docker push nrlabcruk/invar2:${version}

