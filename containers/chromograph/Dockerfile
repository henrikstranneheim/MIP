################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="11"
LABEL software="chromograph"
LABEL software.version="1.1"
LABEL extra.binaries="chromograph"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install pip python=3.7 matplotlib

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

## Download release
RUN wget --no-verbose https://github.com/mikaell/chromograph/archive/refs/tags/v1.1.zip -O /opt/conda/share/chromograph-1.1.zip

## Extract
RUN unzip -o /opt/conda/share/chromograph-1.1.zip -d /opt/conda/share/ && \
    rm /opt/conda/share/chromograph-1.1.zip

## Move to chromograph directory
RUN cd /opt/conda/share/chromograph-1.1 && \
    python -m pip install --no-cache-dir .

WORKDIR /data/
