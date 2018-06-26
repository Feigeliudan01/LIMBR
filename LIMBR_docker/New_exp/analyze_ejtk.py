from LIMBR import simulations

sims = {}

for i in range(1,21):
    analysis = simulations.analyze('standard_' + str(i) + '_true_classes.txt')
    analysis.add_data('standard_' + str(i) + '_LIMBR_processed__jtkout_GammaP.txt','LIMBR')
    analysis.add_data('standard_' + str(i) + '_old_processed__jtkout_GammaP.txt','traditional')
    analysis.add_data('standard_' + str(i) + '_baseline__jtkout_GammaP.txt','baseline')
    analysis.add_data('standard_eigenMS_' + str(i) + '__jtkout_GammaP.txt','eigenMS')
    analysis.calculate_auc()
    sims[i] = analysis.roc_auc