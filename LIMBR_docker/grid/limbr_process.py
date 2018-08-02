from LIMBR import simulations, imputation, batch_fx, old_fashioned

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('twenty_miss_5_NN_'+str(i)+'_with_noise.txt',missingness=0.1,neighbors=5)
    #Impute and Write Output
    to_impute.impute_data('twenty_miss_5_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='twenty_miss_5_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('twenty_miss_5_NN_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('twenty_miss_10_NN_'+str(i)+'_with_noise.txt',missingness=0.1,neighbors=10)
    #Impute and Write Output
    to_impute.impute_data('twenty_miss_10_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='twenty_miss_10_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('twenty_miss_10_NN_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('twenty_miss_15_NN_'+str(i)+'_with_noise.txt',missingness=0.1,neighbors=15)
    #Impute and Write Output
    to_impute.impute_data('twenty_miss_15_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='twenty_miss_15_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('twenty_miss_15_NN_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('thirty_miss_5_NN_'+str(i)+'_with_noise.txt',missingness=0.2,neighbors=5)
    #Impute and Write Output
    to_impute.impute_data('thirty_miss_5_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='thirty_miss_5_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('thirty_miss_5_NN_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('thirty_miss_10_NN_'+str(i)+'_with_noise.txt',missingness=0.2,neighbors=10)
    #Impute and Write Output
    to_impute.impute_data('thirty_miss_10_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='thirty_miss_10_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('thirty_miss_10_NN_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('thirty_miss_15_NN_'+str(i)+'_with_noise.txt',missingness=0.2,neighbors=15)
    #Impute and Write Output
    to_impute.impute_data('thirty_miss_15_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='thirty_miss_15_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('thirty_miss_15_NN_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('forty_miss_5_NN_'+str(i)+'_with_noise.txt',missingness=0.3,neighbors=5)
    #Impute and Write Output
    to_impute.impute_data('forty_miss_5_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='forty_miss_5_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('forty_miss_5_NN_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('forty_miss_10_NN_'+str(i)+'_with_noise.txt',missingness=0.3,neighbors=10)
    #Impute and Write Output
    to_impute.impute_data('forty_miss_10_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='forty_miss_10_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('forty_miss_10_NN_'+str(i)+'_LIMBR_processed.txt')

for i in range(1, 21):
    #Read Raw Data
    to_impute = imputation.imputable('forty_miss_15_NN_'+str(i)+'_with_noise.txt',missingness=0.3,neighbors=15)
    #Impute and Write Output
    to_impute.impute_data('forty_miss_15_NN_'+str(i)+'_imputed.txt')

for i in range(1, 21):
    #Read Imputed Data
    to_sva = batch_fx.sva(filename='forty_miss_15_NN_'+str(i)+'_imputed.txt',design='c',data_type='p',pool='pool_map.p')
    #preprocess data
    to_sva.preprocess_default()
    #perform permutation testing
    to_sva.perm_test(nperm=100)
    #write_output
    to_sva.output_default('forty_miss_15_NN_'+str(i)+'_LIMBR_processed.txt')

