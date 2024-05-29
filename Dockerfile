FROM ubuntu:20.04


ENV TZ="America/Chicago"
ENV LOCATION="America/Chicago"
ENV LC_NUMERIC=en_US.UTF-8
ENV LANGUAGE=en_US
ENV LANG=en_US.UTF-8
ENV LC_ALL=C
RUN apt-get update
RUN apt-get -y install --no-install-recommends git subversion gcc g++ make wget gfortran patch pkg-config file
RUN apt-get -y install --no-install-recommends libgfortran-10-dev libblas-dev liblapack-dev libmetis-dev libnauty2-dev
RUN apt-get -y install --no-install-recommends ca-certificates
RUN apt-get -y install less vim htop
RUN apt-get -y install libopenblas-dev
RUN git clone https://github.com/coin-or/coinbrew /var/coin-or
WORKDIR /var/coin-or
RUN ./coinbrew fetch Alps@master
RUN ./coinbrew --tests none build Alps 
