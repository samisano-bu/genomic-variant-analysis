FROM mcr.microsoft.com/devcontainers/python:3.11

# Install R
RUN apt-get update && \
    apt-get install -y r-base

# Install R packages
RUN R -e "install.packages(c('tidyverse','viridisLite'), repos='https://cloud.r-project.org/')"
