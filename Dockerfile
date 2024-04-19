FROM ubuntu:20.04


ENV TZ="America/Chicago"
RUN apt-get update
RUN apt-get -y install --no-install-recommends git subversion gcc g++ make wget gfortran patch pkg-config file
RUN apt-get -y install --no-install-recommends libgfortran-10-dev libblas-dev liblapack-dev libmetis-dev libnauty2-dev
RUN apt-get -y install --no-install-recommends ca-certificates
RUN apt-get -y install openmpi-bin libopenmpi-dev 
RUN apt-get -y install less vim htop
RUN git clone https://github.com/coin-or/coinbrew /var/coin-or
WORKDIR /var/coin-or
RUN ./coinbrew fetch Alps@master
RUN ./coinbrew build --tests none Alps --enable-static --disable-shared --with-mpi-cflags="$(pkg-config --cflags mpi)" --with-mpi-lflags="$(pkg-config --libs mpi)" MPICC=mpicc MPICXX=mpiCC
RUN apt-get -y install libopenblas-dev
