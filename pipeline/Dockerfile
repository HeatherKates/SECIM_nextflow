# Use the base R image for R version 4.2
FROM rocker/r-ver:4.2.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev libcurl4-openssl-dev libssl-dev pandoc git && \
    apt-get clean

# Set up Bioconductor version compatible with R/4.2
RUN Rscript -e "install.packages('BiocManager')" \
 && Rscript -e "BiocManager::install(version = '3.15')"

# Copy pipeline scripts and resources into the container
WORKDIR /pipeline
COPY . /pipeline

# Install R dependencies from libraries.txt
RUN Rscript -e "pkgs <- scan('/pipeline/libraries.txt', what = 'character', sep = ','); \
    BiocManager::install(pkgs, force = TRUE, ask = FALSE)"

# Copy custom R packages to the Docker image
COPY R-packages/metid /usr/local/lib/R/site-library/metid
COPY R-packages/downloadthis /usr/local/lib/R/site-library/downloadthis
COPY R-packages/MetaboDiff /usr/local/lib/R/site-library/MetaboDiff

# Set the working directory for the container
WORKDIR /workspace

# Specify entry point for R
ENTRYPOINT ["Rscript"]
