from LIMBR import simulations
import pandas as pd

sims = {}

for i in range(1,21):
    analysis = simulations.analyze('standard_' + str(i) + '_true_classes.txt')
    analysis.add_data('standard_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.add_data('standard_' + str(i) + '_old_processed__jtkout_GammaP.txt','traditional')
    analysis.add_data('standard_' + str(i) + '_baseline__jtkout_GammaP.txt','baseline')
    analysis.add_data('standard_eigenMS_' + str(i) + '__jtkout_GammaP.txt','eigenMS')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

data = pd.DataFrame(sims).T
data['Noise'] = 'standard'

sims = {}

for i in range(1,21):
    analysis = simulations.analyze('double_noise_' + str(i) + '_true_classes.txt')
    analysis.add_data('double_noise_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.add_data('double_noise_' + str(i) + '_old_processed__jtkout_GammaP.txt','traditional')
    analysis.add_data('double_noise_' + str(i) + '_baseline__jtkout_GammaP.txt','baseline')
    analysis.add_data('double_noise_eigenMS_' + str(i) + '__jtkout_GammaP.txt','eigenMS')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Noise'] = 'double'
data = pd.concat([data, temp_data])

sims = {}

for i in range(1,21):
    analysis = simulations.analyze('half_noise_' + str(i) + '_true_classes.txt')
    analysis.add_data('half_noise_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.add_data('half_noise_' + str(i) + '_old_processed__jtkout_GammaP.txt','traditional')
    analysis.add_data('half_noise_' + str(i) + '_baseline__jtkout_GammaP.txt','baseline')
    analysis.add_data('half_noise_eigenMS_' + str(i) + '__jtkout_GammaP.txt','eigenMS')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Noise'] = 'half'
data = pd.concat([data, temp_data])

sims = {}

for i in range(1,21):
    analysis = simulations.analyze('tenth_noise_' + str(i) + '_true_classes.txt')
    analysis.add_data('tenth_noise_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.add_data('tenth_noise_' + str(i) + '_old_processed__jtkout_GammaP.txt','traditional')
    analysis.add_data('tenth_noise_' + str(i) + '_baseline__jtkout_GammaP.txt','baseline')
    analysis.add_data('tenth_noise_eigenMS_' + str(i) + '__jtkout_GammaP.txt','eigenMS')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc

temp_data = pd.DataFrame(sims).T
temp_data['Noise'] = 'tenth'
data = pd.concat([data, temp_data])

data.to_csv('aucs.txt',sep='\t')