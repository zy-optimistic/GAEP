# Use Micromamba as a lightweight, fast Conda-enabled base image
FROM mambaorg/micromamba:1.5.1

LABEL maintainer="Bioinformatics User"
LABEL description="Docker image for GAEP: A genome assembly assessment pipeline"

# Set up working directory
WORKDIR /opt/GAEP

# Install system dependencies required for some bioinformatics tools
USER root
RUN apt-get update && apt-get install -y \
    git \
    wget \
    build-essential \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Switch back to micromamba user for environment setup
USER $MAMBA_USER

# Create the environment and install GAEP dependencies
# GAEP relies on several heavy-duty tools: BUSCO, QUAST, BWA, SAMtools, etc.
RUN micromamba create -y -n gaep_env -c conda-forge -c bioconda \
    python=3.8 \
    busco=5.4.7 \
    quast=5.2.0 \
    bwa=0.7.17 \
    samtools=1.17 \
    seqkit=2.4.0 \
    merqury=1.3 \
    meryl=1.3 \
    numpy \
    matplotlib \
    pandas

# Clone the GAEP repository
RUN git clone https://github.com/zy-optimistic/GAEP.git /opt/GAEP

# Set environment variables to ensure the conda environment is active
ENV PATH /opt/conda/envs/gaep_env/bin:$PATH
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Make the GAEP main script executable
RUN chmod +x /opt/GAEP/gaep

# Add the source directory to PATH so GAEP.py can be called directly
ENV PATH /opt/GAEP:$PATH

# Default command
CMD ["python", "/opt/GAEP/gaep", "--help"]
