FROM opencadc/astropy:3.9-slim

RUN apt-get update --no-install-recommends && \
    apt-get dist-upgrade -y && \
    apt-get install -y \
        build-essential \
        git \
        libcfitsio-bin \
        libhdf5-dev \
        wget && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

WORKDIR /usr/src/app

RUN cd /tmp && \
    wget https://support.hdfgroup.org/ftp/HDF5/tools/h5check/src/h5check-2.0.1.tar.gz  && \
    tar xvf h5check-2.0.1.tar.gz && \
    cd h5check-2.0.1 && \
    ./configure && \
    make && \
    cp tool/h5check /usr/local/bin && \
    cd /usr/src/app && \
    rm -rf /tmp/h5check-2.0.1

ARG CAOM2_BRANCH=master
ARG CAOM2_REPO=opencadc
ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc
ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN git clone https://github.com/${CAOM2_REPO}/caom2tools.git && \
    cd caom2tools && \
    git checkout ${CAOM2_BRANCH} && \
    pip install ./caom2utils && \
    cd ..

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${PIPE_REPO}/taosii2caom2@${PIPE_BRANCH}#egg=taosii2caom2

ENTRYPOINT ["/docker-entrypoint.sh"]
