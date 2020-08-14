FROM openanalytics/r-base

LABEL maintainer "Jason Sckrabulis <jason.sckrabulis@gmail.com>"

#system libraries
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
    
#basic shiny libraries
RUN R -e "install.packages(c('shiny','rmarkdown'),repos='https://cloud.r-project.org/')"

#app libraries
Run R -e "install.packages(c(),repos='https://cloud.r-project.org/')"

#copy app to docker image
RUN mkdir /root/euler
COPY euler /root/euler

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/euler')"]
