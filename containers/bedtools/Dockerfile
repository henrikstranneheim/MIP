################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="bedtools"
LABEL software.version="2.30.0"
LABEL extra.binaries="bedtools"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bedtools=2.30.0

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
