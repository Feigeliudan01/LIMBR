from LIMBR import simulations, imputation, batch_fx, old_fashioned

simulation = simulations.simulate()
simulation.generate_pool_map()
simulation.write_output()

#Read Data
to_impute = imputation.imputable('simulated_data_with_noise.txt',0.3)

#Remove Duplicate Peptides
to_impute.deduplicate()

#Drop Rows Over Missing Value Threshold
to_impute.drop_missing()

#Impute and Write Output
to_impute.impute('imputed.txt')



#import data
to_sva = batch_fx.sva(filename='imputed.txt',design='c',data_type='p',pool='pool_map.p')
#normalize for pooled controls
to_sva.pool_normalize()
#calculate timepoints from header
to_sva.get_tpoints()
#calculate correlation with primary trend of interest
to_sva.prim_cor()
#reduce data based on primary trend correlation
to_sva.reduce()
#calculate residuals
to_sva.set_res()
#calculate tks
to_sva.set_tks()
#perform permutation testing
to_sva.perm_test(nperm=100)
#perform eigen trend regression
to_sva.eig_reg(alpha=0.05)
#perform subset svd
to_sva.subset_svd()
#write_output
to_sva.normalize('LIMBR_processed.txt')

#
to_old = old_fashioned.old_fashioned('simulated_data_with_noise.txt',exp_type,pool_map)

#
to_old.pool_normalize()

#
to_old.normalize('old_processed.txt')

#remove unique replicate identifiers to match eJTK header requirements
sed -e 's/_[[:digit:]]//g' LIMBR_processed.txt > temp.txt
cut -f 1,5- temp.txt > LIMBR_processed.txt

sed -e 's/_[[:digit:]]//g' old_processed.txt > temp.txt
cut -f 1,5- temp.txt > old_processed.txt

python2 LIMBR_docker/src/eJTK-CalcP.py -f LIMBR_processed.txt -w LIMBR_docker/src/ref_files/waveform_cosine.txt -a LIMBR_docker/src/ref_files/asymmetries_02-22_by2.txt -s LIMBR_docker/src/ref_files/phases_00-22_by2.txt -p LIMBR_docker/src/ref_files/period24.txt

python2 LIMBR_docker/src/eJTK-CalcP.py -f old_processed.txt -w LIMBR_docker/src/ref_files/waveform_cosine.txt -a LIMBR_docker/src/ref_files/asymmetries_02-22_by2.txt -s LIMBR_docker/src/ref_files/phases_00-22_by2.txt -p LIMBR_docker/src/ref_files/period24.txt