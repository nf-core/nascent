FROM continuumio/miniconda3:4.6.14
MAINTAINER Edmund Miller <Edmund.A.Miller@gmail.com>

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install --yes wget samtools r-essentials bioconductor-deseq2 bioconductor-edger
RUN apt-get update -y
RUN apt-get install -y gcc g++ unzip wget zip
RUN mkdir homer/
RUN curl -fsSL http://homer.ucsd.edu/homer/configureHomer.pl > /homer/configureHomer.pl
RUN perl /homer/configureHomer.pl -install
ENV PATH="/homer/bin/:${PATH}"
RUN perl /homer/configureHomer.pl -install hg18
RUN perl /homer/configureHomer.pl -install hg19
