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
from statsmodels.nonparametric.smoothers_lowess import lowess
from pymc3 import DiscreteUniform, Beta, switch
import pymc3 as pm
import pymc as pm

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
    for i in range(1,len(df.columns)):
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


def prim_cor(arr,l,per,reps):
    cors = []
    for row in arr:
        ave = [np.mean(i) for i in [row[j:j + reps] for j in range(0, len(row), reps)]]
        cors.append((autocorr(ave,per) - autocorr(ave,(per//2))))
    return cors

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

def subset_svd(arr,plist,lam):
    trends = []
    for entry in plist:
        sub = []
        thresh = est_pi_sig(entry,lam)
        if thresh == 'nan':
            return trends
        for i in range(len(entry)):
            if entry[i] < thresh:
                sub.append(arr[i])
        U, s, V = np.linalg.svd(sub)
        trends.append(V.T[:,0])
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

csp_norm = data[data.columns.values[:75]]
wt_norm = data[data.columns.values[75:]]

wt_norm.columns = [i.replace('WT-','CT') for i in wt_norm.columns.values]

wt_norm.columns = [i.replace('-','_') for i in wt_norm.columns.values]

csp_norm.columns = [i.replace('C-','CT') for i in csp_norm.columns.values]

csp_norm.columns = [i.replace('-','_') for i in csp_norm.columns.values]

tpoints = get_tpoints(wt_norm.columns.values)

wt_norm = pd.DataFrame(scale(wt_norm.values,axis=1),columns=wt_norm.columns,index=wt_norm.index)

cors = prim_cor(wt_norm.values,tpoints,12,3)

np.mean(cors)
uncor = [(i<(np.percentile(cors,25))) for i in cors]
len(wt_norm[uncor])
len(wt_norm)
wt_norm_reduced = wt_norm[uncor]

res = get_res(wt_norm_reduced.values,tpoints)

sns.clustermap(res[:1000],col_cluster=False)


tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,tpoints,100,tks_)
end = time.time()
print(end - start)

sigs

ps = eig_reg(wt_norm_reduced.values,res,sigs,.05)

ts = subset_svd(wt_norm_reduced.values,ps,0.5)

plt.plot(range(len(ts[0])),ts[0])

plt.plot(range(len(ts[1])),ts[1])

plt.plot(range(len(ts[2])),ts[2])
plt.plot(range(len(ts[3])),ts[3])

plt.plot(range(len(ts[4])),ts[4])

U, s, V = np.linalg.svd(res)
sig = V.T[:,:len([i for i in itertools.takewhile(lambda x: x < .05, sigs)])]

fin_res = np.dot(linalg.lstsq(np.asarray(sig),wt_norm.values.T)[0].T,np.asarray(sig).T)
fin_res = np.dot(linalg.lstsq(np.asarray(ts).T,wt_norm.values.T)[0].T,np.asarray(ts))

sns.clustermap(fin_res[:100],col_cluster=False)

svd_norm = wt_norm.values - fin_res

svd_norm = pd.DataFrame(svd_norm,index=wt_norm.index,columns=wt_norm.columns)
svd_norm = pd.DataFrame(scale(svd_norm.values,axis=1),columns=svd_norm.columns,index=svd_norm.index)

svd_norm = svd_norm.groupby(level='Protein').mean()


svd_norm.index.names = ['#']

svd_norm.to_csv('wt_reduced_subset_lowess_no_qnorm.txt',sep='\t')

sns.clustermap(svd_norm[:100],col_cluster=False)


csp_norm = pd.DataFrame(scale(csp_norm.values,axis=1),columns=csp_norm.columns,index=csp_norm.index)

tpoints = get_tpoints(wt_norm.columns.values)


cors = prim_cor(csp_norm.values,tpoints,12,3)

sns.distplot(cors)

np.mean(cors)
uncor = [(i<(np.percentile(cors,25))) for i in cors]
len(csp_norm[uncor])
len(csp_norm)
csp_norm_reduced = csp_norm[uncor]

res = get_res(csp_norm_reduced.values,tpoints)

sns.clustermap(csp_norm_reduced[:1000],col_cluster=False)


tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,tpoints,100,tks_)
end = time.time()
print(end - start)

sigs

U, s, V = np.linalg.svd(res)
sig = V.T[:,:len([i for i in itertools.takewhile(lambda x: x < .05, sigs)])]

fin_res = np.dot(linalg.lstsq(np.asarray(sig),wt_norm.values.T)[0].T,np.asarray(sig).T)

sns.clustermap(fin_res[:1000],col_cluster=False)


svd_norm = csp_norm.values - fin_res

svd_norm = pd.DataFrame(svd_norm,index=csp_norm.index,columns=csp_norm.columns)
svd_norm = pd.DataFrame(scale(svd_norm.values,axis=1),columns=svd_norm.columns,index=svd_norm.index)

svd_norm = svd_norm.groupby(level='Protein').mean()


svd_norm.index.names = ['#']

svd_norm.to_csv('csp_lowess.txt',sep='\t')


ps = eig_reg(csp_norm_reduced.values,res,sigs,.05)

ts = subset_svd(wt_norm_reduced.values,ps,0.5)

fin_res = np.dot(linalg.lstsq(np.asarray(ts).T,csp_norm.values.T)[0].T,np.asarray(ts))

sns.clustermap(fin_res[:1000],col_cluster=False)



svd_norm = csp_norm.values - fin_res

svd_norm = pd.DataFrame(svd_norm,index=csp_norm.index,columns=csp_norm.columns)
svd_norm = pd.DataFrame(scale(svd_norm.values,axis=1),columns=svd_norm.columns,index=svd_norm.index)

svd_norm = svd_norm.groupby(level='Protein').mean()


svd_norm.index.names = ['#']

svd_norm.to_csv('csp_reduced_subset_lowess_no_qnorm.txt',sep='\t')

#############################################





np.sum(tks_[:5])
np.sum(tks_[:13])

data = ps[7]

plt.hist(data)



class RandomScanDiscreteMetropolis(pm.step_methods.arraystep.ArrayStep):
    def __init__(self, var, model=None, values=[0,1]):
        model = pm.modelcontext(model)
        self.values = values
        super(RandomScanDiscreteMetropolis, self).__init__([var], [model.fastlogp])

    def astep(self, q0, logp):
        i = np.random.choice(len(q0))
        q = np.copy(q0)
        q[i] = np.random.choice(self.values)
        q_new = pm.step_methods.arraystep.metrop_select(logp(q) - logp(q0), q, q0)
        return q_new


class SequentialScanDiscreteMetropolis(pm.step_methods.arraystep.ArrayStep):
    def __init__(self, var, model=None, values=[0,1]):
        model = pm.modelcontext(model)
        self.values = values
        self.i = 0
        super(SequentialScanDiscreteMetropolis, self).__init__([var], [model.fastlogp])

    def astep(self, q0, logp):
        q = np.copy(q0)
        q[self.i] = np.random.choice(self.values)
        self.i = (self.i + 1) % len(q)
        q_new = pm.step_methods.arraystep.metrop_select(logp(q) - logp(q0), q, q0)
        return q_new

p = pm.Uniform('p',0.,1.)

assignment = pm.Categorical('assignment',[p,1-p],size=len(data))


mu = pm.Beta('mu',[.5,.5],[1.,1.],size=2)
etas = pm.Uniform('etas', 0, 1, size=2)

@pm.deterministic
def a_i(assignment=assignment,mu=mu,etas=etas):
    return mu[assignment]*etas[assignment]

@pm.deterministic
def b_i(assignment=assignment,mu=mu,etas=etas):
    return mu[assignment]*(1-etas[assignment])

observations = pm.Beta('obs',alpha=a_i,beta=b_i,value=data,observed=True)

model = pm.Model([p,assignment,a_i,b_i])

mcmc = pm.MCMC(model)

mcmc.sample(100000,1000,10)

a = mcmc.trace('a_i')[:]
b = mcmc.trace('b_i')[:]
ass = mcmc.trace('assignment')[:]
p = mcmc.trace('p')[:]
len(a)

plt.plot(p[:])

plt.plot(ass[:,1])

plt.hist(b[:,0])
plt.hist(b[:,1])
plt.hist(a[:,1])
plt.hist(a[:,0])

plt.plot(a[:,0])

plt.plot(b[:,0])
plt.plot(b[:,1])

beta = stats.beta
y1 = beta.pdf(np.linspace(0,1,1000),.53,.48)
plt.plot(np.linspace(0,1,1000),y1)

len(a)
np.shape(a)
plt.hist(a)

with pm.Model() as model:
    #priors
    p = pm.Uniform( "p", 0 , 1) #this is the fraction that come from mean1 vs mean2
    ber = pm.Bernoulli( "ber", p = p, shape=len(data)) # produces 1 with proportion p.
    sd = pm.Uniform('sd', lower=0, upper=20)
    means = pm.Normal('means', mu=[0, 0], sd=15, shape=2)
    process = pm.Beta('process',mu=means[ber],sd=sd,observed=data)

with model:
    step1 = pm.Metropolis(vars=[p])
    step2 = pm.Metropolis(vars=[sd, means])
    step3 = pm.BinaryMetropolis([ber], scaling=.01)
    trace = pm.sample(10000, [step1, step2])

with model:
    step1 = pm.Metropolis(vars=[p])
    step2 = pm.Metropolis(vars=[sd, means])
    step3 = SequentialScanDiscreteMetropolis(var=ber, values=[0,1])
    tr = pm.sample(10000, step=[step1, step2] + [step3]*len(data))

pm.traceplot(trace[1000:])
plt.show()

import pymc as pm
etas = pm.Normal('sigmas', mu=0.1, tau=1000, size=2)
centers = pm.Normal('centers', [0.3, 0.7], [1/(0.1)**2, 1/(0.1)**2], size=2)
alpha  = pm.Beta('alpha', alpha=2, beta=3)

@pm.deterministic
def a_i(assignment=assignment,mu=mu,sigmas=sigmas):
    return mu[assignment]*sigmas[assignment]

@pm.deterministic
def b_i(assignment=assignment,mu=mu,sigmas=sigmas):
    return mu[assignment]*(1-sigmas[assignment])

category = pm.Container([pm.Categorical("category%i" % i, [alpha, 1 - alpha]) for i in range(nsamples)])
observations = pm.Container([pm.Normal('samples_model%i' % i, mu=centers[category[i]], tau=1/(sigmas[category[i]]**2), value=samples[i], observed=True) for i in range(nsamples)])
