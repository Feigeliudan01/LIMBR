LIMBR: Learning and Imputation for Mass-spec Bias Reduction
===========================================================

LIMBR provides a streamlined tool set for imputation of missing data followed by modelling and removal of batch effects.  The software was designed for proteomics datasets, with an emphasis on circadian 
proteomics data, but can be applied to any time course or blocked experiments which produce large amounts of data, such as RNAseq. The two main classes are imputable, which performs missing data imputation, and sva, which performs 
modelling and removal of batch effects.

----------
Motivation
----------

Decreasing costs and increasing ambition are resulting in larger Mass-spec (MS) experiments.  MS experiments have a few limitations which are exacerbated by this increasing scale, namely batch effects and missing data.  Many 
downstream statistical analyses require complete cases for analysis, however, MS produces some missing data at random meaning that as the number of experiments increase the number of peptides rejected due to missing data actually 
*increases*.  This is obviously not good, but fortunately there is a solution!  If the missing data for observations missing only a small number of data points are imputed this issue can be overcome and that's the first thing that 
LIMBR does.  The second issue with large scaler MS experiments is batch effects.  As the number of samples increases, the number of batches necessary for sample processing also increases.  Batch effects from sample processing are 
known to have a large effect on MS data and increasing the number of batches means more batch effects and a higher proportion of observations affected by at least one batch effect.  Here LIMBR capitolizes on the larger amount of 
data and the known correlation structure of the data set to model these batch effects so that they can be removed.

--------
Features
--------

* KNN based imputation of missing data.

* SVD based modelling and removal of batch effects.

* Built for circadian and non-circadian time series as well as block designs

-------------
Example Usage
-------------

------------
Installation
------------

pip install limbr

-----------
How to Use?
-----------

----
TODO
----

* Switch to long format files for greater interoperability and more easily specified file format.

* Add unit tests to docstrings where possible.

-------
Credits
-------

The sva based methods build on work for micro-array datasets by Jeffrey Leek, with particular reliance on his PhD Thesis from the University of Washington.  Other contributions in the field include EigenMS.

-------
License
-------

Copyright 2017 Alexander M. Crowell

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

