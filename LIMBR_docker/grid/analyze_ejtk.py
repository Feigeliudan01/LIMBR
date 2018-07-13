from LIMBR import simulations
import pandas as pd

sims = {}

for i in range(1,21):
    analysis = simulations.analyze('twenty_miss_5_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('twenty_miss_5_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

data = pd.DataFrame(sims).T
data['Params'] = 'twenty_miss_5_NN'

sims = {}

for i in range(1,21):
    analysis = simulations.analyze('twenty_miss_10_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('twenty_miss_10_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Params'] = 'twenty_miss_10_NN'
data = pd.concat([data, temp_data])

for i in range(1,21):
    analysis = simulations.analyze('twenty_miss_15_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('twenty_miss_15_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Params'] = 'twenty_miss_15_NN'
data = pd.concat([data, temp_data])

for i in range(1,21):
    analysis = simulations.analyze('thirty_miss_5_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('thirty_miss_5_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Params'] = 'thirty_miss_5_NN'
data = pd.concat([data, temp_data])

for i in range(1,21):
    analysis = simulations.analyze('thirty_miss_10_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('thirty_miss_10_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Params'] = 'thirty_miss_10_NN'
data = pd.concat([data, temp_data])

for i in range(1,21):
    analysis = simulations.analyze('thirty_miss_15_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('thirty_miss_15_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Params'] = 'thirty_miss_15_NN'
data = pd.concat([data, temp_data])

for i in range(1,21):
    analysis = simulations.analyze('forty_miss_5_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('forty_miss_5_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Params'] = 'forty_miss_5_NN'
data = pd.concat([data, temp_data])

for i in range(1,21):
    analysis = simulations.analyze('forty_miss_10_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('forty_miss_10_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Params'] = 'forty_miss_10_NN'
data = pd.concat([data, temp_data])

for i in range(1,21):
    analysis = simulations.analyze('forty_miss_15_NN_' + str(i) + '_true_classes.txt')
    analysis.add_data('forty_miss_15_NN_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Params'] = 'forty_miss_15_NN'
data = pd.concat([data, temp_data])

data.to_csv('aucs.txt',sep='\t')
