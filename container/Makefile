
version ?= 1.0.3

build:
	docker build --tag nrlabcruk/invar2:${version} .

release: build
	docker push nrlabcruk/invar2:${version}

