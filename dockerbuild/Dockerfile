FROM acrowell/limbr:v1
MAINTAINER Alec Crowell <alexander.m.crowell.gr@dartmouth.edu>
 
# Install all the software needed to run the pipeline
RUN apt-get -qq update
RUN sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | sudo apt-key add -
RUN apt-get install -y wget git python python3 nano python-tk python3-pip python-pip python3-tk r-base r-base-dev
RUN su - -c "R -e \"install.packages('circular', repos = 'http://cran.rstudio.com/')\""
RUN pip2 install --upgrade pip
RUN pip3 install --upgrade pip
RUN pip3 install -U pandas scipy sklearn numpy statsmodels cython matplotlib tqdm multiprocess seaborn astroML rpy2
RUN pip2 install -U numpy scipy cython matplotlib pandas statsmodels

# clone the most recent version of the pipeline
RUN git clone https://github.com/aleccrowell/LIMBR.git
RUN mkdir ./LIMBR/data/
RUN mkdir ./LIMBR/output/actual/
RUN cp ./data/*.txt ./LIMBR/data/
RUN mkdir ./LIMBR/output/simdata/
RUN mkdir ./LIMBR/output/figs/
WORKDIR /LIMBR/src/bin/
RUN python3 setup.py build_ext --inplace
RUN python setup.py build_ext --inplace
WORKDIR /LIMBR/

#make actual data files
RUN make -f make_sim
RUN make -f make_actual
RUN make -f make_figs
