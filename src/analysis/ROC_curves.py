import pandas as pd
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt

l = pd.read_csv('output/simdata/circ_lowess_classifications.txt',sep='\t',header=None)

scores = l[0].values.astype(float)
y = l[1].values.astype(int)


fpr, tpr, thresholds = metrics.roc_curve(y, scores)

roc_auc = metrics.auc(fpr, tpr)

plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)

plt.show()
scores = np.asarray([1-i for i in l[0].values.astype(float)])
fpr, tpr, thresholds = metrics.roc_curve(y, scores)
roc_auc = metrics.auc(fpr, tpr)
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
plt.show()

b = pd.read_csv('output/simdata/block_classifications.txt',sep='\t',header=None)
scoresb = np.asarray([1-i for i in b[0].values.astype(float)])
yb = b[1].values.astype(int)
fprb, tprb, thresholdsb = metrics.roc_curve(yb, scoresb)
roc_aucb = metrics.auc(fprb, tprb)
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)

plt.plot(fprb, tprb, label='ROC curve (area = %0.2f)' % roc_aucb)

plt.legend(loc="lower right")

plt.show()
