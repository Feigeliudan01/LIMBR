FROM acrowell/limbr:v1
MAINTAINER Alec Crowell <alexander.m.crowell.gr@dartmouth.edu>
 
# Install all the software needed to run the pipeline
RUN apt-get -qq update
RUN apt-get install -y wget git python nano python-tk python-pip
RUN pip install --upgrade pip
RUN pip install -U pandas scipy sklearn numpy statsmodels cython matplotlib tqdm
  
# clone the most recent version of the pipeline
RUN git clone https://github.com/aleccrowell/LIMBR.git
RUN mkdir ./LIMBR/data/
RUN mkdir ./LIMBR/output/actual/
RUN cp ./data/*.txt ./LIMBR/data/
WORKDIR /LIMBR/src/bin/
RUN python setup.py build_ext --inplace
WORKDIR /LIMBR/

#make actual data files
RUN make -f make_actual

