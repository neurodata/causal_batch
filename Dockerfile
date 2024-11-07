FROM bioconductor/bioconductor_docker:devel

RUN Rscript -e "install.packages(c('tidyverse', 'parallelDist', 'MatchIt', 'energy', 'multcomp', 'survey', 'dplyr'))"
RUN Rscript -e "BiocManager::install('sva')"
RUN Rscript -e "install.packages('igraph')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('mgcv')"
RUN Rscript -e "install.packages('entropy')"
RUN Rscript -e "install.packages(c('mltools'))"
RUN Rscript -e "install.packages(c('cdcsis'))"
RUN apt-get update
RUN apt-get install -y python3-venv
RUN mkdir ~/.virtualenvs
RUN python3 -m venv ~/.virtualenvs/neuroharm
RUN git clone https://github.com/rpomponio/neuroHarmonize.git /neuroharm

RUN . ~/.virtualenvs/neuroharm/bin/activate && pip install neuroCombat==0.2.12 neuroHarmonize numpy
RUN Rscript -e "install.packages('reticulate')"

RUN . ~/.virtualenvs/neuroharm/bin/activate && pip install pandas
RUN . ~/.virtualenvs/neuroharm/bin/activate && pip install statsmodels nibabel
RUN Rscript -e "install.packages(c('GeneralisedCovarianceMeasure', 'weightedGCM'))"
RUN Rscript -e "install.packages(c('momentchi2', 'MASS', 'devtools'))"
RUN Rscript -e "library(devtools); install_github('ericstrobl/RCIT')"
RUN Rscript -e "library(devtools); install_github('neurodata/causal_batch')"
