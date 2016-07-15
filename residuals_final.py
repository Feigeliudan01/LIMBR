%matplotlib inline
from random import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
from numpy.linalg import svd
import time
from sklearn.preprocessing import scale, MinMaxScaler
from sklearn.decomposition import RandomizedPCA
from pylab import *
from math import atan2
from scipy.stats import linregress
import itertools



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

def qnorm(df):
    ref = pd.concat([df[col].sort_values().reset_index(drop=True) for col in df], axis=1, ignore_index=True).mean(axis=1).values
    for i in range(0,len(df.columns)):
        df = df.sort_values(df.columns[i])
        df[df.columns[i]] = ref
    return df.sort_index()

def get_res(arr,l):
    res = []
    for row in arr:
        ys = lowess(row, l,delta=4)[:,1]
        res.append(row - ys)
    return res

def get_tpoints(l):
    tpoints = [i.replace('CT','') for i in l]
    tpoints = [int(i.split('_')[0]) for i in tpoints]
    return np.asarray(tpoints)

def autocorr(l,shift):
    return dot(l, np.roll(l, shift)) / dot(l, l)


def prim_cor(arr,l,per):
    cors = []
    for row in arr:
        ave = []
        for k in set(l):
            ave.append((np.mean([row[i] for i, j in enumerate(l) if j == k])*1000000))
        cors.append((autocorr(ave,per) - autocorr(ave,(per//2))))
    return np.asarray(cors)

def eig_reg(arr,resarr,perm,a):
    U, s, V = np.linalg.svd(resarr)
    sig = V.T[:,:len([i for i in itertools.takewhile(lambda x: x < a, perm)])]
    pvals = []
    for trend in sig.T:
        temp = []
        for row in arr:
            slope, intercept, r_value, p_value, std_err = linregress(row,trend)
            temp.append(p_value)
        pvals.append(temp)
    return pvals

def normalize(arr,resarr,perm,a):
    pca = RandomizedPCA(n_components=len([i for i in perm if i <= a]))
    pca.fit(resarr)
    return arr  - pca.inverse_transform(pca.transform(resarr))

def est_pi_naught(probs_naught,lam):
    return len([i for i in probs_naught if i > lam])/(len(probs_naught)*(1-lam))

def est_pi_sig(probs_sig,l):
    pi_0 = est_pi_naught(probs_sig,l)
    if pi_0 > 1:
        return 'nan'
    sp = np.sort(probs_sig)
    return sp[int(floor((1-pi_0)*len(probs_sig)))]

def subset_svd(arr,plist,lam,resarr):
    _, _, bt = np.linalg.svd(resarr)
    trends = []
    for j, entry in enumerate(plist):
        sub = []
        thresh = est_pi_sig(entry,lam)
        if thresh == 'nan':
            return trends
        for i in range(len(entry)):
            if entry[i] < thresh:
                sub.append(arr[i])
        U, s, V = np.linalg.svd(sub)
        temp = []
        for trend in V:
            _, _, _, p_value, _ = linregress(bt[j],trend)
            temp.append(p_value)
        trends.append(V.T[:,np.argmin(temp)])
    return trends

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

def get_tks(resarr):
    pca = RandomizedPCA()
    pca.fit(resarr)
    return pca.explained_variance_ratio_

os.chdir('/Users/user4574/Dropbox (Dunlap-Loros labs)/MSLearn')

data = pd.read_csv('imputed_peptide_final.txt',sep='\t')
data = data.set_index(['Peptide','Protein'])
norm_map = gen_norm_dict(data.columns.values)
data = pool_normalize(data,norm_map)
data = data.replace([np.inf, -np.inf], np.nan)
data = data.dropna()
data = data.sort_index(axis=1)
data = qnorm(data)

len(wt_norm)

csp_norm = data[data.columns.values[:75]]
wt_norm = data[data.columns.values[75:]]

wt_norm.columns = [i.replace('WT-','CT') for i in wt_norm.columns.values]

wt_norm.columns = [i.replace('-','_') for i in wt_norm.columns.values]

wt_norm = wt_norm[wt_norm.columns.values[3:]]

csp_norm.columns = [i.replace('C-','CT') for i in csp_norm.columns.values]

csp_norm.columns = [i.replace('-','_') for i in csp_norm.columns.values]

csp_norm = csp_norm[csp_norm.columns.values[3:]]

wt_norm = pd.DataFrame(scale(wt_norm.values,axis=1),columns=wt_norm.columns,index=wt_norm.index)

csp_norm = pd.DataFrame(scale(csp_norm.values,axis=1),columns=csp_norm.columns,index=csp_norm.index)

tpoints = get_tpoints(wt_norm.columns.values)

cors = prim_cor(wt_norm.values,tpoints,12)

uncor = [(i<(np.percentile(cors,25))) for i in cors]

wt_norm_reduced = wt_norm[uncor]

res = get_res(wt_norm_reduced.values,tpoints)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,tpoints,100,tks_)
end = time.time()
print(end - start)

sigs

ps = eig_reg(wt_norm_reduced.values,res,sigs,.05)

ts = subset_svd(wt_norm_reduced.values,ps,0.5, res)

fin_res = np.dot(linalg.lstsq(np.asarray(ts).T,wt_norm.values.T)[0].T,np.asarray(ts))

sns.clustermap(fin_res[:1000],col_cluster=False)

svd_norm = wt_norm.values - fin_res

svd_norm = pd.DataFrame(svd_norm,index=wt_norm.index,columns=wt_norm.columns)
svd_norm = pd.DataFrame(scale(svd_norm.values,axis=1),columns=svd_norm.columns,index=svd_norm.index)

svd_norm = svd_norm.groupby(level='Protein').mean()


svd_norm.index.names = ['#']

svd_norm.to_csv('wt_reduced_subset_lowess_qnorm_final_cor.txt',sep='\t')


fin_res = pd.DataFrame(fin_res,index=wt_norm.index,columns=wt_norm.columns)
fin_res = fin_res.groupby(level='Protein').mean()
fin_res.index.names = ['#']
fin_res.to_csv('wt_residuals_cor.txt',sep='\t')



tpoints = get_tpoints(csp_norm.columns.values)


cors = prim_cor(csp_norm.values,tpoints,12)

uncor = [(i<(np.percentile(cors,25))) for i in cors]

csp_norm_reduced = csp_norm[uncor]

res = get_res(csp_norm_reduced.values,tpoints)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,tpoints,100,tks_)
end = time.time()
print(end - start)

sigs

ps = eig_reg(csp_norm_reduced.values,res,sigs,.05)

ts = subset_svd(csp_norm_reduced.values,ps,0.5, res)

fin_res = np.dot(linalg.lstsq(np.asarray(ts).T,csp_norm.values.T)[0].T,np.asarray(ts))

svd_norm = csp_norm.values - fin_res

svd_norm = pd.DataFrame(svd_norm,index=csp_norm.index,columns=csp_norm.columns)
svd_norm = pd.DataFrame(scale(svd_norm.values,axis=1),columns=svd_norm.columns,index=svd_norm.index)

svd_norm = svd_norm.groupby(level='Protein').mean()


svd_norm.index.names = ['#']

svd_norm.to_csv('csp_reduced_subset_lowess_qnorm_final_cor.txt',sep='\t')

sns.clustermap(fin_res[:1000],col_cluster=False)

fin_res = pd.DataFrame(fin_res,index=csp_norm.index,columns=csp_norm.columns)
fin_res = fin_res.groupby(level='Protein').mean()
fin_res.index.names = ['#']
fin_res.to_csv('csp_residuals_cor.txt',sep='\t')


##############


data = pd.read_csv('Jen_rnaseq_Alec_raw_counts.csv')


data = data.set_index('Unnamed: 0')
data.index.names = ['#']


data.columns = [i.replace('set','') for i in data.columns.values]
data.columns = ['CT' + i.split("-")[1] + '_' + i.split("-")[0] for i in data.columns.values]

data = data.sort_index(axis=1)
data = qnorm(data)

data = pd.DataFrame(scale(data.values,axis=1),columns=data.columns,index=data.index)

tpoints = get_tpoints(data.columns.values)

cors = prim_cor(data.values,tpoints,12)

data = data.drop(data.index[[i for i, j in enumerate(cors) if np.isnan(j) == True]])

cors = prim_cor(data.values,tpoints,12)

uncor = [(i<(np.percentile(cors,25))) for i in cors]

data_reduced = data[uncor]

res = get_res(data_reduced.values,tpoints)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,tpoints,100,tks_)
end = time.time()
print(end - start)

sigs

ps = eig_reg(data_reduced.values,res,sigs,.05)

ts = subset_svd(data_reduced.values,ps,0.5, res)

fin_res = np.dot(linalg.lstsq(np.asarray(ts).T,data.values.T)[0].T,np.asarray(ts))

svd_norm = data.values - fin_res

svd_norm = pd.DataFrame(svd_norm,index=data.index,columns=data.columns)
svd_norm = pd.DataFrame(scale(svd_norm.values,axis=1),columns=svd_norm.columns,index=svd_norm.index)

svd_norm.index.names = ['#']

svd_norm.to_csv('rna_reduced_subset_lowess_qnorm_final_cor.txt',sep='\t')
