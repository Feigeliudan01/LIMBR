FROM ubuntu:latest
MAINTAINER Alec Crowell <alexander.m.crowell@gmail.com>
#builds an interactive ubuntu environment in docker for running r2r simulations

# Install all the software needed to run the pipeline
ENV DEBIAN_FRONTEND noninteractive
ENV LANG C.UTF-8
RUN apt-get -qq update
RUN apt-get -y install software-properties-common curl
RUN apt-get install -y wget git python python3 nano python-tk python3-pip python-pip python3-tk
RUN pip2 install --upgrade pip
RUN pip3 install --upgrade pip
RUN pip3 install -U numpy cython matplotlib seaborn astroML jinja2 limbr
RUN pip2 install -U numpy cython matplotlib limbr

RUN git clone https://github.com/aleccrowell/LIMBR.git
RUN mv LIMBR temp
RUN mv temp/LIMBR_docker LIMBR
RUN rm -rf temp
WORKDIR /LIMBR/src/bin/
RUN python3 setup.py build_ext --inplace
RUN python setup.py build_ext --inplace

WORKDIR /LIMBR/New_exp
RUN python3 Generate_data.py
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' standard_$i\_baseline.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do cut -f 1,5- temp_$i.txt > standard_$i\_baseline.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f New_exp/standard_$i\_baseline.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/New_exp
RUN rm temp_*.txt
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' double_noise_$i\_baseline.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do cut -f 1,5- temp_$i.txt > double_noise_$i\_baseline.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f New_exp/double_noise_$i\_baseline.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/New_exp
RUN rm temp_*.txt
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' half_noise_$i\_baseline.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do cut -f 1,5- temp_$i.txt > half_noise_$i\_baseline.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f New_exp/half_noise_$i\_baseline.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/New_exp
RUN rm temp_*.txt
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' tenth_noise_$i\_baseline.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do cut -f 1,5- temp_$i.txt > tenth_noise_$i\_baseline.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f New_exp/tenth_noise_$i\_baseline.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/New_exp
RUN python3 limbr_process.py

