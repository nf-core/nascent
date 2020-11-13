SHELL:=/usr/bin/env bash
NXF_VER:=20.10.0

resume: install
	./nextflow -d run main.nf -profile test,docker -resume

clean_test: install
	./nextflow -d run main.nf -profile test,docker

install: ./nextflow

./nextflow:
	export NXF_VER="$(NXF_VER)" && \
	curl -fsSL get.nextflow.io | bash

clean:
	rm -f .nextflow.log*
	rm -rf .nextflow*
	rm -rf work


run: install
	./nextflow run main.nf

