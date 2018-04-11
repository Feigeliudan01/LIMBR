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
LIMBR does.  The second issue with larger scale MS experiments is batch effects.  As the number of samples increases, the number of batches necessary for sample processing also increases.  Batch effects from sample processing are 
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

```
from LIMBR import simulations, imputation, batch_fx

simulation = simulations.simulate()
simulation.generate_pool_map()
simulation.write_output()

#Read Data
to_impute = imputable.imputable('simulated_data_with_noise.txt',0.3)
#Remove Duplicate Peptides
to_impute.deduplicate()
#Drop Rows Over Missing Value Threshold
to_impute.drop_missing()
#Impute and Write Output
to_impute.impute('imputed.txt')

#set SVA parameters
#raw data
data_path = 'imputed.txt'
#circadian experimental design
exp_design = 'c'
#Proteomic data type
exp_type = 'p'
#Sample blocks (Used in non timecourse designs)
blocks = None
#Pool Normalization map
pool_map = 'pool_map.p'
#percentage data reduction
psub = 25
#number of permutations
nperm = 10
#number of processors
num_proc = 1
#alpha (significance level)
a= 0.05
#lambda (subset threshold)
l = 0.5

#import data
to_sva = sva.sva(data_path,exp_design,exp_type,blocks,pool_map)
#normalize for pooled controls
to_sva.pool_normalize()
#calculate timepoints from header
to_sva.get_tpoints()
#calculate correlation with primary trend of interest
to_sva.prim_cor()
#reduce data based on primary trend correlation
to_sva.reduce(psub)
#calculate residuals
to_sva.set_res()
#calculate tks
to_sva.set_tks()
#perform permutation testing
to_sva.perm_test(nperm,num_proc)
#perform eigen trend regression
to_sva.eig_reg(a)
#perform subset svd
to_sva.subset_svd(l)
#write_output
to_sva.normalize('LIMBR_processed.txt')
```

------------
Installation
------------

pip install limbr

-------------
API Reference
-------------

http://limbr.readthedocs.io/en/latest/

-----------
How to Use?
-----------

----
TODO
----

* Switch to long format files for greater interoperability and more easily specified file format.

* Add unit tests to docstrings where possible.

* Review ensuring maximum Vectorization

-------
Credits
-------

K nearest neighbors as an imputation method was originally proposed by Gustavo Batista in 2002 (http://conteudo.icmc.usp.br/pessoas/gbatista/files/his2002.pdf) and has seen a great deal of success since.

The sva based methods build on work for micro-array datasets by Jeffrey Leek, with particular reliance on his PhD Thesis from the University of Washington (https://digital.lib.washington.edu/researchworks/bitstream/handle/1773/9586/3290558.pdf?sequence=1).

-------
License
-------

© 2017 Alexander M. Crowell: BSD-3
