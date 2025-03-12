FROM neurodata/causal_batch:0.1.0

RUN apt-get update
RUN apt-get install -y python3
RUN apt-get install -y g++ gcc
RUN apt-get install -y python3-venv

RUN pip3 install jupyter jupyterlab
RUN python3 -m venv /root/.virtualenvs/causal

RUN /root/.virtualenvs/causal/bin/pip install \
    seaborn dask dask[distributed] \
    git+https://github.com/py-why/dodiscover@6e64456 \
    statsmodels rpy2

RUN R -e "install.packages('GeneralisedCovarianceMeasure', dependencies=TRUE)"
RUN R -e "install.packages('weightedGCM', dependencies=TRUE)"
RUN apt-get -y install libcurl4-openssl-dev

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev

RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('energy')"
RUN R -e "install.packages('vegan')"
RUN R -e "install.packages('cdcsis')"
RUN R -e "install.packages('reticulate')"


RUN R -e "library(devtools); install_github('ebridge2/RCIT')"
RUN R -e "library(devtools); install_github('neurodata/causal_batch@7d927b4')"
RUN /root/.virtualenvs/causal/bin/pip install git+https://github.com/ebridge2/hyppo@a8fe83f
