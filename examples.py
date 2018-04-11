from LIMBR import simulations, imputable, sva

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
to_impute.impute('sva_imputed.txt')

#set SVA parameters

#raw data
data_path = 'sva_imputed.txt'

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
nperm = 1000

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