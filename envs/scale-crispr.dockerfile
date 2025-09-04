FROM continuumio/miniconda3
# Path to yaml file containing conda dependencies
ARG CONDA_YML_PATH=.
RUN apt-get clean && apt-get update && \
 DEBIAN_FRONTEND=noninteractive apt-get -y install \
 build-essential \
 curl libz-dev \
 unzip && \
rm -rf /var/lib/apt/lists/*

# Install the conda environment
COPY $CONDA_YML_PATH/scale-crispr.conda.yml /environment.yml
RUN conda env create --quiet -f /environment.yml && \
 conda clean -a --yes
# Instead of 'conda activate'
ENV PATH="/opt/conda/envs/scale-crispr/bin:${PATH}"

# bcParser, STAR, etc.
COPY $CONDA_YML_PATH/download-scale-tools.sh /
RUN /download-scale-tools.sh /tools
ENV PATH="/tools:${PATH}"
