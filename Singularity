From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Ignacio Tripodi, Margaret Gruca
    DESCRIPTION Singularity image containing all requirements for the nf-core/nascent pipeline
    VERSION 1.0

%environment
    PATH=/opt/conda/envs/nf-core-nascent-1.0dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
