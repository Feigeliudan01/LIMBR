from LIMBR import simulations, imputation, batch_fx, old_fashioned

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('standard_'+str(i)+'_with_noise.txt',0.3)
    #Impute and Write Output
    to_impute.impute_data('standard_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='standard_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('standard_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    to_old = old_fashioned.old_fashioned(filename='standard_'+str(i)+'_with_noise.txt',data_type='p',pool='pool_map.p')
    to_old.pool_normalize()
    to_old.normalize('standard_'+str(i)+'_old_processed.txt')