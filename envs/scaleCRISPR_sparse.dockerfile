FROM continuumio/miniconda3
# Path to yaml file containing conda dependencies
ARG CONDA_YML_PATH=.

# Install the conda environment
COPY $CONDA_YML_PATH/scaleCRISPR_sparse.conda.yml /environment.yml
RUN conda env create --quiet -f /environment.yml && \
 conda clean -a --yes
# Instead of 'conda activate'
ENV PATH="/opt/conda/envs/scaleCRISPR_sparse/bin:${PATH}"