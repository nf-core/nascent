FROM nfcore/base
LABEL authors="Ignacio Tripodi (ignacio.tripodi@colorado.edu), Margaret Gruca (margaret.gruca@colorado.edu)" \
      description="Docker image containing all requirements for nf-core/nascent pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-nascent-1.0/bin:$PATH
