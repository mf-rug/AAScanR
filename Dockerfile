# AAScanR docker file
# https://github.com/mf-rug/AAScanR

# get shiny server and R from the rocker project
FROM rocker/shiny:4.2.1

# use amd64 platform as I'm building this on Mac M1, check if this affects when run
FROM --platform=linux/amd64 ubuntu

# system libraries for httr
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev
  

# install R packages required for AAScanR
RUN R -e 'install.packages(c(\
              "shiny", \
              "shinyjs", \
              "DT", \
              "ggplot2", \
              "magrittr", \
              "stringr", \
              "readr", \
              "dplyr", \
              "tidyr", \
              "tools", \
              "reticulate", \
              "zip" \
            ), \
            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2023-02-01"\
          )'
          
          # copy the app directory into the image
          
COPY ./app/* /srv/shiny-server/

# run app
CMD ["/usr/bin/shiny-server"]
