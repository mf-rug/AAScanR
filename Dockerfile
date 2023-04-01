# AAScanR docker file
# https://github.com/mf-rug/AAScanR

# get shiny server and R from the rocker project
FROM rocker/shiny:4.2.2

# use amd64 platform as I'm building this on Mac M1, check if this affects when run
#FROM --platform=linux/amd64 ubuntu

# system libraries for httr
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev
    
# python for reticulate
RUN apt-get install -y \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv
    
RUN apt update 
RUN apt-get -y install git 

RUN pip3 install Bio warnings random pyramid numpy scipy

RUN git clone 'https://github.com/matteoferla/DirEvo_tools'
RUN cd DirEvo_tools
RUN python3 setup.py install
    
# install R packages required for AAScanR
RUN ls -la /usr/local/bin

RUN /usr/local/bin/R -e 'install.packages(c(\
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
          
#RUN /usr/local/bin/R -e 'reticulate::use_python("/usr/bin/python3, required = T)")'
#RUN /usr/local/bin/R -e 'reticulate::install_miniconda()'

COPY ./app/* /srv/shiny-server/

#start shiny-server as the shiny user rather than root
USER shiny

ENV SHINY_LOG_STDERR=1

# run app
CMD ["/usr/bin/shiny-server"]
