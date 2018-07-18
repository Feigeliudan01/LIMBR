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

WORKDIR /LIMBR/grid
RUN python3 Generate_data.py
RUN python3 limbr_process.py

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' twenty_miss_5_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt twenty_miss_5_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/twenty_miss_5_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' twenty_miss_10_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt twenty_miss_10_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/twenty_miss_10_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' twenty_miss_15_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt twenty_miss_15_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/twenty_miss_15_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' thirty_miss_5_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt thirty_miss_5_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/thirty_miss_5_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' thirty_miss_10_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt thirty_miss_10_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/thirty_miss_10_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' thirty_miss_15_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt thirty_miss_15_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/thirty_miss_15_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' forty_miss_5_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt forty_miss_5_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/forty_miss_5_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' forty_miss_10_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt forty_miss_10_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/forty_miss_10_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

WORKDIR /LIMBR/grid
RUN for i in $(seq 1 20); do sed -e 's/_[[:digit:]]//g' forty_miss_15_NN_$i\_LIMBR_processed.txt  > temp_$i.txt; done
RUN for i in $(seq 1 20); do mv temp_$i.txt forty_miss_15_NN_$i\_LIMBR_processed.txt; done
WORKDIR /LIMBR/
RUN for i in  $(seq 1 20); do python2 src/eJTK-CalcP.py -f grid/forty_miss_15_NN_$i\_LIMBR_processed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

COPY ./analyze_ejtk.py /LIMBR/grid/analyze_ejtk.py
RUN pip install -U limbr
WORKDIR /LIMBR/grid
RUN python3 analyze_ejtk.py
