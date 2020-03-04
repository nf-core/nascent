FROM nfcore/base:1.9
LABEL authors="Ignacio Tripodi, Margaret Gruca" \
      description="Docker image containing all software requirements for the nf-core/nascent pipeline"

# Set the locale (needed for IGVtools and some Python libraries)
#RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
#    locale-gen
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8  

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-nascent-1.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-nascent-1.1 > nf-core-nascent-1.1.yml
