FROM rocker/tidyverse:4.1.0

# Install ViteRbi

RUN R -e "devtools::install_github('tf2/ViteRbi', upgrade = 'never')"

# Install linux packages required for karyoploteR

RUN apt update
RUN apt install -y \
  # for Rhtslib
  libbz2-dev liblzma-dev

# Install karyoploteR

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

RUN R -e "BiocManager::install('karyoploteR')"
