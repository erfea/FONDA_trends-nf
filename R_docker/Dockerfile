###########################################################################
# :AUTHOR: Katarzyna Ewa Lewinska  [lewinska@hu-berlin.de]
# :DATE: 14 March 2023
#
# :Description: A Docker image for a container to run remotePARTS
#
# :Updates:
#
# :2Do:
#
# :Disclaimer:  The author of this code accepts no responsibility for errors or omissions in this work
#               and shall not be liable for any damage caused by these.
###########################################################################

FROM r-base:4.2.2

# disable interactive frontends
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
  apt-get install -y \
  libcurl4-openssl-dev \
  libfontconfig1-dev \
  libfreetype6-dev \
  libfribidi-dev \
  libgdal-dev \
  libharfbuzz-dev \
  libjpeg-dev \
  libpng-dev \
  libssh2-1-dev \
  libssl-dev \
  libtiff5-dev \
  libxml2-dev \
  zlib1g-dev \
  cmake
RUN Rscript -e "install.packages(c('devtools', 'dplyr'))" && \
    Rscript -e "devtools::install_github('morrowcj/remotePARTS')" && \
    Rscript -e "install.packages(c('snow', 'parallel', 'doParallel', 'foreach', 'graphics'))" && \
    Rscript -e "install.packages(c('raster', 'data.table', 'rgdal'))" && \

# RUN R -e "install.packages(c('devtools','testthat', 'roxygen2'), dependencies=TRUE, repos='https://cran.rstudio.com')" && \
# RUN R -e 'install.packages("devtools", dependencies=TRUE, repos="http://cran.rstudio.com/")'  && \
#    R -e ', install.packages("remotePARTS", dependencies=TRUE, repos="http://cran.rstudio.com/")'  && \

#RUN Rscript -e install.packages("devtools")  && \
#    Rscript -e install.packages("remotePARTS")  && \

    # Clear installation data
    apt-get clean && rm -r /var/cache/

ENV HOME /home/docker
# ENV PATH $HOME/docker/bin

WORKDIR /home/docker

CMD ['Rsc']
