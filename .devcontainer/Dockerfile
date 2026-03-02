FROM mcr.microsoft.com/devcontainers/python:3.11

# System deps for common R packages (readr/stringr/tidyr/ggplot2, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install only the R packages you actually use
RUN R -e "install.packages(c('dplyr','ggplot2','readr','tidyr','stringr','forcats','viridisLite','scales'), repos='https://cloud.r-project.org/', Ncpus=parallel::detectCores())"
