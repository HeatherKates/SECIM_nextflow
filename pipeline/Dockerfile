FROM rocker/r-ver:4.3.1

RUN apt-get update && apt-get install -y \
    libxml2-dev libcurl4-openssl-dev libssl-dev pandoc git && \
    apt-get clean

# Install MetaboAnalystR dependencies
RUN Rscript -e "install.packages(c('devtools', 'rmarkdown', 'BiocManager'))"
RUN Rscript -e "BiocManager::install(c('MetaboAnalystR'))"

# Copy custom scripts (if needed)
WORKDIR /pipeline
COPY R/ ./R/

ENTRYPOINT ["Rscript"]

