FROM ubuntu:latest
MAINTAINER Alec Crowell <alexander.m.crowell@gmail.com>
#builds an interactive ubuntu environment in docker for running r2r simulations

# Install all the software needed to run the pipeline
ENV DEBIAN_FRONTEND noninteractive
ENV LANG C.UTF-8
RUN apt-get -qq update
RUN apt-get -y install software-properties-common curl
RUN apt-get install -y wget git python python3 nano python-tk python3-pip python-pip python3-tk
RUN add-apt-repository ppa:git-core/ppa
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install git-lfs
RUN git lfs install
RUN pip2 install --upgrade pip
RUN pip3 install --upgrade pip
RUN pip3 install -U numpy cython matplotlib seaborn astroML jinja2 limbr
RUN pip2 install -U numpy cython matplotlib limbr

RUN git lfs clone https://github.com/aleccrowell/LIMBR.git
RUN mv LIMBR temp
RUN mv temp/LIMBR_docker LIMBR
RUN rm -rf temp
WORKDIR /LIMBR/src/bin/
RUN python3 setup.py build_ext --inplace
RUN python setup.py build_ext --inplace

WORKDIR /LIMBR/New_exp
RUN python3 Generate_data.py
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do sed -e 's/_[[:digit:]]//g' standard_$i_baseline.txt  > temp_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do cut -f 1,5- temp_$i.txt > standard_$i_baseline.txt; done
RUN rm temp_*.txt
WORKDIR /LIMBR/
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python2 src/eJTK-CalcP.py -f New_exp/standard_$i_baseline.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done