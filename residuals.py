import pandas as pd
import numpy as np
import os
from numpy.linalg import svd
from fbpca import pca

os.chdir('/Users/user4574/Dropbox (Dunlap-Loros labs)/MSLearn')

data = pd.read_csv('imputed_peptide_final.txt',sep='\t')
data = data.set_index(['Peptide','Protein'])
data.head()

#need to take list of classes iterate through the list and for the indexes sharing a class average each row and subtract that row from those columns

#then need to do svd on the resulting residuals matrix do significance permutation tests etc...


def get_res(arr,l):
    m = {}
    for v in set(l):
        indices = [i for i, x in enumerate(l) if x == v]
        m[v] = np.mean(arr[:,indices],axis=1)
    ma = np.zeros(np.shape(arr))
    for i in range(len(l)):
        ma[:,i]=m[l[i]]
    return np.subtract(arr,ma)

def get_tks(resarr):
     U, s, V = pca(resarr, raw=True)
     denom = np.sum([j * j for j in s])
     tk = []
     for i in range(len(s)):
         tk.append((s[i]*s[i])/denom)
     return tk

def perm_test(resarr,l,n,tks):
    rstar = resarr
    out = np.zeros(len(tks))
    for j in range(n):
        for i in range(rstar.shape[0]):
            np.random.shuffle(rstar[i,:])
        resstar = get_res(rstar,l)
        tkstar = get_tks(resstar)
        for m in range(len(tks)):
            if tkstar[m] >= tks[m]:
                out[m] += 1
    return out/n

len(data.columns)
classes = np.random.randint(2, size=len(data.columns))
res = get_res(data.values,classes)
classes
data.values

tks = get_tks(res)
np.sum(tks)
perm_test(res,classes,10,tks)
