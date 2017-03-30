FROM acrowell/limbr:v1
MAINTAINER Alec Crowell <alexander.m.crowell.gr@dartmouth.edu>
 
# Install all the software needed to run the pipeline
RUN apt-get -qq update
RUN sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -
RUN apt-get install -y wget git python python3 nano python-tk python3-pip python-pip python3-tk r-base r-base-dev
RUN su - -c "R -e \"install.packages('circular', repos = 'http://cran.rstudio.com/')\""
RUN pip2 install --upgrade pip
RUN pip3 install --upgrade pip
RUN pip3 install -U pandas scipy sklearn numpy statsmodels cython matplotlib tqdm multiprocess seaborn astroML rpy2
RUN pip2 install -U numpy scipy cython matplotlib pandas statsmodels
ENV TLDIR texlive
RUN mkdir -p $TLDIR
WORKDIR $TLDIR

# Download and uncompress latest version of texlive
RUN wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz && \
   tar -xzf install-tl-unx.tar.gz

# Run the installer
RUN cd install-tl-*/. && touch default.profile && \
   ./install-tl --persistent-downloads --profile default.profile

# Set symlinks to installed program
RUN ln -s /usr/local/texlive/2015/bin/* /opt/texbin

# Set the PATH to use the newly installed tex
ENV PATH $PATH:/opt/texbin

# Clean installation files
WORKDIR /
RUN rm -rf $TLDIR


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
