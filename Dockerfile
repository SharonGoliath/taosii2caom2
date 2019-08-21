FROM python:3.6-alpine

RUN apk --no-cache add \
        bash \
        coreutils \
        gcc \
        git \
        g++ \
        libffi-dev \
        libmagic \
        libxml2-dev \
        libxslt-dev \
        make \
        musl-dev \
        openssl-dev

RUN pip install aenum && \
        pip install astropy && \
        pip install cadcdata && \
        pip install cadctap && \
        pip install caom2repo && \
        pip install funcsigs && \
        pip install future && \
        pip install numpy && \
        pip install PyYAML && \
        pip install spherical-geometry && \
        pip install vos && \
        pip install xml-compare

WORKDIR /usr/src/app
RUN git clone https://github.com/opencadc-metadata-curation/caom2tools.git && \
  cd caom2tools && git pull origin master && \
  pip install ./caom2utils && pip install ./caom2pipe && cd ..

# the make check fails :( twice now, so I know it's not a random
# error

#  ./configure --prefix=/usr/local/hdf5  && \
# the hdf5.h file ends up not where pip is looking :(
RUN wget https://s3.amazonaws.com/hdf-wordpress-1/wp-content/uploads/manual/HDF5/HDF5_1_10_5/source/hdf5-1.10.5.tar.gz && \
  tar -xzvf hdf5-1.10.5.tar.gz && cd hdf5-1.10.5 && \
  ./configure  && \
  make && make check && make install && make check-install

RUN pip install --no-binary=h5py h5py

RUN git clone https://github.com/opencadc-metadata-curation/taosii2caom2.git && \
  pip install ./taosii2caom2

COPY ./docker-entrypoint.sh /

ENTRYPOINT ["/docker-entrypoint.sh"]
