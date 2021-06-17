FROM rocker/tidyverse:4.1.0

# Install ViteRbi
RUN R -e "devtools::install_github('tf2/ViteRbi', upgrade = 'never')"
