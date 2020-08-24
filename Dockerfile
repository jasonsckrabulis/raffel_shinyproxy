FROM openanalytics/r-base

LABEL maintainer "Jason Sckrabulis <jason.sckrabulis@gmail.com>"

# system dependencies
RUN apt-get update && apt-get install -y \
   sudo \
   pandoc \
   pandoc-citeproc \
   libcurl4-gnutls-dev \
   libcairo2-dev \
   libxt-dev \
   libssl-dev \
   libssh2-1-dev \
   libssl1.0.0

# shiny libraries
RUN R -e "install.packages(c('shiny','rmarkdown'), repos='https://cloud.r-project.org/')"

# app libraries
RUN R -e "install.packages(c('ggplot2','cowplot','gridExtra','grid','knitr'),repos="https://cloud.r-project.org/')"

# copy app to image
RUN mkdir /root/TPCapp
COPY TPCapp /root/TPCapp

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/TPCapp')"]
