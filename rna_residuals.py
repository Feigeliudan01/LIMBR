%matplotlib inline
import pandas as pd
import numpy as np
import os
from numpy.linalg import svd
from fbpca import pca
import time
from sklearn.preprocessing import scale, MinMaxScaler
from sklearn.decomposition import RandomizedPCA
from scipy.optimize import curve_fit
from pylab import *
from math import atan2


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

def get_tpoints(l):
    tpoints = [i.replace('CT','') for i in l]
    tpoints = [int(i.split('_')[0]) for i in tpoints]
    return np.asarray(tpoints)

def get_res(arr,l):
    m = {}
    for v in set(l):
        indices = [i for i, x in enumerate(l) if x == v]
        m[v] = np.mean(arr[:,indices],axis=1)
    ma = np.zeros(np.shape(arr))
    for i in range(len(l)):
        ma[:,i]=m[l[i]]
    return np.subtract(arr,ma)

def fitSine(tList,yList,freq):
   '''
       freq in Hz
       tList in sec
   returns
       phase in degrees
   '''
   b = matrix(yList).T
   rows = [ [sin(freq*2*pi*t), cos(freq*2*pi*t), 1] for t in tList]
   A = matrix(rows)
   (w,residuals,rank,sing_vals) = lstsq(A,b)
   phase = atan2(w[1,0],w[0,0])*180/pi
   amplitude = norm([w[0,0],w[1,0]],2)
   bias = w[2,0]
   return (phase,amplitude,bias)

def get_res(arr,l):
    res = []
    for row in arr:
        def harm(x, p1,p2,p3,p4):
            return p1*np.cos(2*np.pi*p2*x + 2*np.pi*p3) + p4
        amplitude = wt_norm.values[0].max() - wt_norm.values[0].min()
        popt, pcov = curve_fit(harm, tpoints, wt_norm.values[0], p0=(amplitude,.043478261,0,1), bounds=([0,.043478260,-23,0], [3., .043478262, 23,1.5]))
        res.append(row - harm(l,popt[0],popt[1],popt[2],popt[3]))
    return res

def get_res(arr,l,frequency):
    res = []
    for pep in arr:
        (phaseEst,amplitudeEst,biasEst) = fitSine(l,pep,frequency)
        res.append(pep - amplitudeEst*sin(l*frequency*2*pi+phaseEst*pi/180.0)+biasEst)
    return res

#def get_tks(resarr):
#     U, s, V = pca(resarr, raw=True)
     #U, s, V = np.linalg.svd(resarr)
#     denom = np.sum([j * j for j in s])
#     tk = []
#     for i in range(len(s)):
#         tk.append((s[i]*s[i])/denom)
#     return tk

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
data = pool_normalize(data,norm_map)
data = data.replace([np.inf, -np.inf], np.nan)
data = data.dropna()
#data = qnorm(data)
#data = np.log2(np.add(data,1.0))
data = data.sort_index(axis=1)
data = rescale_peptides(data)

data = pd.read_csv('Jen_rnaseq_Alec_raw_counts.csv')
data = data.set_index(['transcript'])



data.columns = ['CT'+i.split('-')[1]+'_'+i.split('-')[0][-1] for i in data.columns]
data = data.sort_index(axis=1)

data.head()

data.columns

classes = [j for i in range(1,14) for j in [i]*3] + [14,14,15,15] + [j for i in range(16,22) for j in [i]*3] + [22,22,23,23,24,24]


csp_norm = data[data.columns.values[:75]]
wt_norm = data[data.columns.values[75:]]

wt_norm.head()

wt_norm.columns = [i.replace('WT-','CT') for i in wt_norm.columns.values]

wt_norm.columns = [i.replace('-','_') for i in wt_norm.columns.values]

wt_norm = wt_norm[wt_norm.columns.values[3:]]

csp_norm.columns = [i.replace('C-','CT') for i in csp_norm.columns.values]

csp_norm.columns = [i.replace('-','_') for i in csp_norm.columns.values]

csp_norm = csp_norm[csp_norm.columns.values[3:]]

(phaseEst,amplitudeEst,biasEst) = fitSine(wt_norm.values[0],tpoints,.043478261)

amplitudeEst

plt.scatter(tpoints,wt_norm.values[0])
plt.plot(tpoints,(amplitudeEst*sin(tpoints*.043478261*2*pi+phaseEst*pi/180.0)+biasEst))

def harm(x, p1,p2,p3,p4):
    return p1*np.cos(2*np.pi*p2*x + 2*np.pi*p3) + p4
amplitude = wt_norm.values[0].max() - wt_norm.values[0].min()
popt, pcov = curve_fit(harm, tpoints, wt_norm.values[0], p0=(amplitude,.043478261,0,1), bounds=([0,.043478260,-23,0], [3., .043478262, 23,1.5]))

amplitude

plt.scatter(tpoints,wt_norm.values[0])
plt.plot(tpoints,harm(tpoints,popt[0],popt[1],popt[2],popt[3]))

data.head()

#exp_class = [0] * 3 + [j for i in range(1,12) for j in [i]*3]*2 + [j for i in range(1,3) for j in [i]*3] + [12] * 20

#classes = exp_class + [i + 13 for i in exp_class]


exp_class = [j for i in range(1,12) for j in [i]*3]*2 + [j for i in range(1,3) for j in [i]*3]

#classes = exp_class + [i + 12 for i in exp_class]

#exp_class = [j for i in range(1,25) for j in [i]*3]

classes = exp_class

len(classes)

len(csp_norm.columns)

classes

#need to take list of classes iterate through the list and for the indexes sharing a class average each row and subtract that row from those columns

#then need to do svd on the resulting residuals matrix do significance permutation tests etc...



tpoints = get_tpoints(csp_norm.columns.values)

res = get_res(csp_norm.values,classes)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,classes,100,tks_)
end = time.time()
print(end - start)

sigs

svd_norm = normalize(csp_norm.values,res,sigs,0.05)


svd_norm = pd.DataFrame(svd_norm,index=csp_norm.index,columns=csp_norm.columns)

svd_norm = svd_norm.groupby(level='Protein').mean()

svd_norm.index.names = ['#']

svd_norm.to_csv('stnd_circ_csp_svd_norm.txt',sep='\t')



tpoints = get_tpoints(wt_norm.columns.values)

res = get_res(wt_norm.values,classes)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,classes,100,tks_)
end = time.time()
print(end - start)

sigs

svd_norm = normalize(wt_norm.values,res,sigs,0.05)


svd_norm = pd.DataFrame(svd_norm,index=wt_norm.index,columns=wt_norm.columns)

svd_norm = svd_norm.groupby(level='Protein').mean()

svd_norm.index.names = ['#']

svd_norm.to_csv('stnd_circ_wt_svd_norm.txt',sep='\t')



svd_norm.index.names = ['#']

svd_norm.to_csv('rep_rna_svd_norm.txt',sep='\t')


svd_norm = normalize(wt_norm.values,res,sigs,0.05)


svd_norm = pd.DataFrame(svd_norm,index=wt_norm.index,columns=wt_norm.columns)

#pnorm = pool_normalize(svd_norm,norm_map)

#pnorm.to_csv('SVD_normalized.csv')
svd_norm = svd_norm.groupby(level='Protein').mean()

svd_norm.index.names = ['#']

svd_norm.to_csv('new_csp_svd_norm.txt',sep='\t')

svd_norm.head()

#svd_norm = 2**svd_norm

svd_norm = svd_norm.groupby(level='Protein').mean()

csp_norm = svd_norm[svd_norm.columns.values[:75]]
wt_norm = svd_norm[svd_norm.columns.values[75:]]

wt_norm.columns = [i.replace('WT-','CT') for i in wt_norm.columns.values]

wt_norm.columns = [i.replace('-','_') for i in wt_norm.columns.values]

wt_norm.index.names = ['#']

wt_norm.to_csv('new_wt_svd_norm.txt',sep='\t')
wt_norm.head()

csp_norm.columns = [i.replace('C-','CT') for i in csp_norm.columns.values]

csp_norm.columns = [i.replace('-','_') for i in csp_norm.columns.values]

csp_norm.index.names = ['#']

csp_norm.to_csv('new_csp_svd_norm.txt',sep='\t')
csp_norm.head()

out.head()
