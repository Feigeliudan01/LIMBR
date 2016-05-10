import pandas as pd
import numpy as np
import os
from numpy.linalg import svd
from fbpca import pca
import time
from sklearn.preprocessing import scale, MinMaxScaler
from sklearn.decomposition import RandomizedPCA
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

pandas2ri.activate()

def gen_norm_dict(l):
    newd = {}
    for i in range(len(l)):
        newd[l[i]] = int(np.ceil((i+1)/10))
    return newd

def pool_normalize(df,dmap):
    newdf = pd.DataFrame(index=df.index)
    for column in df.columns.values:
        if column[0] == 'W' and column[4] != 'o':
            newdf[column] = df[column].div(df['WT_pool_'+'%02d' % dmap[column]],axis='index')
        if column[0] == 'C' and column[4] != 'o':
            newdf[column] = df[column].div(df['C_pool_'+'%02d' % dmap[column]],axis='index')
    return newdf

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
     #U, s, V = np.linalg.svd(resarr)
     denom = np.sum([j * j for j in s])
     tk = []
     for i in range(len(s)):
         tk.append((s[i]*s[i])/denom)
     return tk

def get_tks(resarr):
    pca = RandomizedPCA()
    pca.fit(resarr)
    return pca.explained_variance_ratio_

def perm_test(resarr,l,n,tks):
    rstar = np.copy(resarr)
    out = np.zeros(len(tks))
    for j in range(n):
        for i in range(rstar.shape[0]):
            np.random.shuffle(rstar[i,:])
        resstar = get_res(rstar,l)
        tkstar = get_tks(resstar)
        #tkstar = get_tks(rstar)
        for m in range(len(tks)):
            if tkstar[m] > tks[m]:
                out[m] += 1
    return out/n

def qnorm(df):
    ref = pd.concat([df[col].sort_values().reset_index(drop=True) for col in df], axis=1, ignore_index=True).mean(axis=1).values
    for i in range(1,len(df.columns)):
        df = df.sort_values(df.columns[i])
        df[df.columns[i]] = ref
    return df.sort_index()

def rescale_peptides(df):
    vals = scale(df.values.astype(float),axis=1)
    newdf = pd.DataFrame(vals,index=df.index,columns=df.columns)
    return newdf

def normalize(arr,resarr,perm,a):
    pca = RandomizedPCA(n_components=len([i for i in perm if i <= a]))
    pca.fit(resarr)
    return arr  - pca.inverse_transform(pca.transform(resarr))


os.chdir('/Users/user4574/Dropbox (Dunlap-Loros labs)/MSLearn')

data = pd.read_csv('imputed_peptide_final.txt',sep='\t')
data = data.set_index(['Peptide','Protein'])
norm_map = gen_norm_dict(data.columns.values)
data = qnorm(data)
data = np.log2(np.add(data,1.0))
data = data.sort_index(axis=1)

exp_class = [0] * 3 + [j for i in range(1,12) for j in [i]*3]*2 + [j for i in range(1,3) for j in [i]*3] + [12] * 20

classes = exp_class + [i + 13 for i in exp_class]

#need to take list of classes iterate through the list and for the indexes sharing a class average each row and subtract that row from those columns

#then need to do svd on the resulting residuals matrix do significance permutation tests etc...

res = get_res(data.values,classes)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,classes,20,tks_)
end = time.time()
print(end - start)

sigs


svd_norm = normalize(data.values,res,sigs,0.05)


svd_norm = pd.DataFrame(svd_norm,index=data.index,columns=data.columns)

pnorm = pool_normalize(svd_norm,norm_map)

pnorm.to_csv('SVD_normalized.csv')



2**pnorm

pnorm.groupby(level='Protein').mean().to_csv('SVD_normalized.csv')

pnorm.shape

pnorm.head()

out.head()
