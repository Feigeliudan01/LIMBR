import pandas as pd
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt

l = pd.read_csv('output/simdata/circ_lowess_1_classifications.txt',sep='\t',header=None)

scores = l[0].values.astype(float)
y = l[1].values.astype(int)

scores = np.asarray([1-i for i in l[0].values.astype(float)])
fpr, tpr, thresholds = metrics.roc_curve(y, scores)
roc_auc = metrics.auc(fpr, tpr)

b = pd.read_csv('output/simdata/block_1_classifications.txt',sep='\t',header=None)
scoresb = np.asarray([1-i for i in b[0].values.astype(float)])
yb = b[1].values.astype(int)
fprb, tprb, thresholdsb = metrics.roc_curve(yb, scoresb)
roc_aucb = metrics.auc(fprb, tprb)

bl = pd.read_csv('output/simdata/simulated_data_baseline_1_classifications.txt',sep='\t',header=None)
scoresbl = np.asarray([1-i for i in bl[0].values.astype(float)])
ybl = bl[1].values.astype(int)
fprbl, tprbl, thresholdsbl = metrics.roc_curve(ybl, scoresbl)
roc_aucbl = metrics.auc(fprbl, tprbl)

n = pd.read_csv('output/simdata/simulated_data_with_noise_1_classifications.txt',sep='\t',header=None)
scoresn = np.asarray([1-i for i in n[0].values.astype(float)])
yn = n[1].values.astype(int)
fprn, tprn, thresholdsn = metrics.roc_curve(yn, scoresn)
roc_aucn = metrics.auc(fprn, tprn)

plt.plot(fpr, tpr, label=' Lowess ROC curve (area = %0.2f)' % roc_auc, color = 'y')

plt.plot(fprb, tprb, label='Block ROC curve (area = %0.2f)' % roc_aucb, color='b')

plt.plot(fprbl, tprbl, label='Baseline ROC curve (area = %0.2f)' % roc_aucbl, color='g')

plt.plot(fprn, tprn, label='Noise ROC curve (area = %0.2f)' % roc_aucn, color='r')

plt.legend(loc="lower right")
plt.axis([0, 1, 0, 1])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.savefig('output/simdata/ROC_curves_1.pdf')

l = pd.read_csv('output/simdata/circ_lowess_2_classifications.txt',sep='\t',header=None)

scores = l[0].values.astype(float)
y = l[1].values.astype(int)

scores = np.asarray([1-i for i in l[0].values.astype(float)])
fpr, tpr, thresholds = metrics.roc_curve(y, scores)
roc_auc = metrics.auc(fpr, tpr)

b = pd.read_csv('output/simdata/block_2_classifications.txt',sep='\t',header=None)
scoresb = np.asarray([1-i for i in b[0].values.astype(float)])
yb = b[1].values.astype(int)
fprb, tprb, thresholdsb = metrics.roc_curve(yb, scoresb)
roc_aucb = metrics.auc(fprb, tprb)

bl = pd.read_csv('output/simdata/simulated_data_baseline_2_classifications.txt',sep='\t',header=None)
scoresbl = np.asarray([1-i for i in bl[0].values.astype(float)])
ybl = bl[1].values.astype(int)
fprbl, tprbl, thresholdsbl = metrics.roc_curve(ybl, scoresbl)
roc_aucbl = metrics.auc(fprbl, tprbl)

n = pd.read_csv('output/simdata/simulated_data_with_noise_2_classifications.txt',sep='\t',header=None)
scoresn = np.asarray([1-i for i in n[0].values.astype(float)])
yn = n[1].values.astype(int)
fprn, tprn, thresholdsn = metrics.roc_curve(yn, scoresn)
roc_aucn = metrics.auc(fprn, tprn)

plt.plot(fpr, tpr, label=' Lowess ROC curve (area = %0.2f)' % roc_auc, color = 'y')

plt.plot(fprb, tprb, label='Block ROC curve (area = %0.2f)' % roc_aucb, color='b')

plt.plot(fprbl, tprbl, label='Baseline ROC curve (area = %0.2f)' % roc_aucbl, color='g')

plt.plot(fprn, tprn, label='Noise ROC curve (area = %0.2f)' % roc_aucn, color='r')

plt.legend(loc="lower right")
plt.axis([0, 1, 0, 1])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.savefig('output/simdata/ROC_curves_2.pdf')

l = pd.read_csv('output/simdata/circ_lowess_3_classifications.txt',sep='\t',header=None)

scores = l[0].values.astype(float)
y = l[1].values.astype(int)

scores = np.asarray([1-i for i in l[0].values.astype(float)])
fpr, tpr, thresholds = metrics.roc_curve(y, scores)
roc_auc = metrics.auc(fpr, tpr)

b = pd.read_csv('output/simdata/block_3_classifications.txt',sep='\t',header=None)
scoresb = np.asarray([1-i for i in b[0].values.astype(float)])
yb = b[1].values.astype(int)
fprb, tprb, thresholdsb = metrics.roc_curve(yb, scoresb)
roc_aucb = metrics.auc(fprb, tprb)

bl = pd.read_csv('output/simdata/simulated_data_baseline_3_classifications.txt',sep='\t',header=None)
scoresbl = np.asarray([1-i for i in bl[0].values.astype(float)])
ybl = bl[1].values.astype(int)
fprbl, tprbl, thresholdsbl = metrics.roc_curve(ybl, scoresbl)
roc_aucbl = metrics.auc(fprbl, tprbl)

n = pd.read_csv('output/simdata/simulated_data_with_noise_3_classifications.txt',sep='\t',header=None)
scoresn = np.asarray([1-i for i in n[0].values.astype(float)])
yn = n[1].values.astype(int)
fprn, tprn, thresholdsn = metrics.roc_curve(yn, scoresn)
roc_aucn = metrics.auc(fprn, tprn)

plt.plot(fpr, tpr, label=' Lowess ROC curve (area = %0.2f)' % roc_auc, color = 'y')

plt.plot(fprb, tprb, label='Block ROC curve (area = %0.2f)' % roc_aucb, color='b')

plt.plot(fprbl, tprbl, label='Baseline ROC curve (area = %0.2f)' % roc_aucbl, color='g')

plt.plot(fprn, tprn, label='Noise ROC curve (area = %0.2f)' % roc_aucn, color='r')

plt.legend(loc="lower right")
plt.axis([0, 1, 0, 1])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.savefig('output/simdata/ROC_curves_3.pdf')
