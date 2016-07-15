%matplotlib inline
from __future__ import division
from matplotlib import pyplot as plt
import numpy as np
import pymc3 as pm
import scipy as sp
import seaborn as sns
from statsmodels.datasets import get_rdataset
from theano import tensor as T
from operator import mul
from pymc3 import *
from functools import reduce
import math as m
import pandas as pd
#theano.config.compute_test_value = 'warn'
theano.config.compute_test_value = 'off'
import theano
import time

brna = pd.read_csv('rna_reduced_subset_lowess_qnorm_final_cor__jtkout_GammaP.txt',sep='\t')
brna = brna.set_index('ID')

bcprot = pd.read_csv('csp_reduced_subset_lowess_qnorm_final_cor__jtkout_GammaP.txt',sep='\t')
bcprot = bcprot.set_index('ID')


bwprot = pd.read_csv('wt_reduced_subset_lowess_qnorm_final_cor__jtkout_GammaP.txt',sep='\t')
bwprot = bwprot.set_index('ID')

x = np.reshape(bwprot['GammaP'].values,(1,-1))
y = np.reshape(bcprot['GammaP'].values,(1,-1))



with pm.Model() as model:
    a = pm.Gamma('a',10.,10.,shape=4)
    def dens(x,y):
        N = 1000
        m = T.dscalar('m')
        n = T.dscalar('n')
        linmax = lambda m,n: T.maximum((m+n-1),0)
        low = linmax(x,y)
        k = T.dscalar('k')
        l = T.dscalar('l')
        linmin = lambda k,l: T.minimum(k,l)
        up = linmin(x,y)
        u = T.dvector('u')
        myint = lambda u: ((u**(a[3]-1))*((x-u)**(a[2]-1))*((y-u)**(a[1]-1))*((1-x-y+u)**(a[0]-1)))
        num = T.dscalar('num')
        ma = T.dscalar('ma')
        mi = T.dscalar('mi')
        tlinspace = lambda ma,mi,num: T.dot(((T.reshape(T.arange(num),(-1,1)) + .5) / num),(ma-mi))+T.extra_ops.repeat(mi,num,axis=0)
        lin = tlinspace(up,low,N)
        area = T.sum(myint(lin))*(up-low)/N
        return T.log(T.gamma(T.sum(a))/T.prod(T.gamma(a))*area)
    Bibeta = pm.DensityDist("Bibeta", dens,  observed={"x": x, 'y':y})

with model:
    start = {}
    start["a"] = [80,10,10,1]
    step = pm.Metropolis(vars=[a])
    trace = pm.sample(2000, step=step, start= start)


with model:
    step = pm.Metropolis(vars=[a])
    trace = pm.sample(2000, step=step)

with model:
    start = pm.find_MAP()

start

with model:
    start = pm.find_MAP()
    trace = pm.sample(2000, start=start)

pm.traceplot(trace, varnames=["a"])

np.mean(trace['a'][::50,0])
N= 1000
delta = 0.05
xgrid = np.arange(delta, 1.0, delta)
ygrid = np.arange(delta, 1.0, delta)
X, Y = np.meshgrid(xgrid, ygrid)
np.shape(X)
a = [84.193916952007285,0.16078286199088215,0.20358242471646312,0.1276771579030774]
a = [60.,20.,20.,1.]
Z = []

start = time.time()
for i in xgrid:
    temp = []
    for j in ygrid:
        temp.append(calc_dens(a,i,j,N))
    Z.append(temp)
end = time.time()
print(end - start)

#need to write a theano function that estimates the DensityDist
#good news is that NUTS with find_MAP seems to converge instantly

def calc_dens(a, x, y, N):
    m = T.dscalar('m')
    n = T.dscalar('n')
    linmax = theano.function([m,n],T.maximum((m+n-1),0))
    low = linmax(x,y)
    k = T.dscalar('k')
    l = T.dscalar('l')
    linmin = theano.function([k,l],T.minimum(k,l))
    up = linmin(x,y)
    u = T.dvector('u')
    myint = theano.function([u],((u**(a[3]-1))*((x-u)**(a[2]-1))*((y-u)**(a[1]-1))*((1-x-y+u)**(a[0]-1))))
    num = T.dscalar('num')
    ma = T.dscalar('ma')
    mi = T.dscalar('mi')
    #tlinspace = theano.function([ma,mi,num],T.dot(((T.reshape(T.arange(num),(-1,1)) + .5) / num),(ma-mi))+T.extra_ops.repeat(mi,num,axis=0))
    tlinspace = theano.function([ma,mi,num],((T.arange(num) + .5) / num)*(ma-mi)+mi)
    lin = tlinspace(up,low,N)
    intout = myint(lin)
    mix = T.dvector('mix')
    auc = T.dscalar('auc')
    quad =T.dvector('quad')
    upb = T.dscalar('upb')
    lowb = T.dscalar('lowb')
    calc = theano.function([mix, quad, upb,lowb,num], T.log(T.gamma(T.sum(mix))/T.prod(T.gamma(mix))*T.sum(quad)*(upb-lowb)/num))
    return calc(a,intout,up,low,N)


start = time.time()
calc_dens(a,.5,.5,1000)
end = time.time()
print(end - start)

imgplot = plt.imshow(Z, vmin=np.min(Z), vmax=np.max(Z), origin='lower', extent=[0,1,0,1])
imgplot.set_cmap('spectral')

sns.jointplot(x=np.reshape(bwprot['GammaP'].values[:1000],(1,-1)), y = np.reshape(brna['GammaP'].values[:1000],(1,-1)), kind='kde')
sns.jointplot(x=x, y = y)




with pm.Model() as model:
    p = pm.Uniform( "p", 0 , 1)
    component = pm.Bernoulli( "component", p = p, shape=N) # produces 1 with proportion p.
    mean = pm.Uniform('mean', 0., 1., shape=K)
    non_null_mu = pm.Deterministic('mu', mean[component])
    mu = pm.switch(T.eq(0, component), .5, non_null_mu)
    sigma = pm.Uniform('sigma', 0., 0.5, shape=K)
    non_null_sd = pm.Deterministic('sd', sigma[component])
    sd = pm.switch(T.eq(0, component), np.sqrt(1/12), non_null_sd)
    obs = pm.Beta('obs', mu=mu, sd=sd, observed=data)





with model:
    start = {}
    start["component"] = [0]*500 + [1]*500
    start["p"] = 0.5
    start["sigma"] = np.array([.1,.2])
    start["mean"] = np.array([.5, .5])

    step1 = pm.Slice([p, sigma, mean])
    step2 = pm.BinaryMetropolis([component], scaling=.01)
    trace = pm.sample(40000, start=start, step=[step1, step2])
    pm.traceplot(trace[20000:][::50], varnames=["p","sigma","mean"])
    print([int(model.logp(t)) for t in trace[20000:][::400]])
    print([trace[i]['component'].sum() for i in np.arange(0,len(trace),400)])
x_plot = np.linspace(0, 1, 500)

fig, ax = plt.subplots(figsize=(8, 6))

def calc_alph_bet(m,s):
    a =  -1*(m*((s**2) + (m ** 2)-m))/(s ** 2)
    b = (((s**2)+(m ** 2) - m)*(m-1))/(s ** 2)
    return a, b


np.mean(trace['mean'][40000:][:,0])
np.mean(trace['sigma'][40000:][:,0])

np.mean(trace['mean'][20000:][:,1])
np.mean(trace['sigma'][20000:][:,1])

calc_alph_bet(1/2,.288666)
calc_alph_bet(0.57655532865125192,0.2435294992952412)
calc_alph_bet(0.21195605432652997,0.22794638804124415)

n_bins = 20
ax.hist(old_faithful_df.std_waiting, bins=n_bins, color=blue, lw=0, alpha=0.5);

ax.set_xlabel('Standardized waiting time between eruptions');
ax.set_ylabel('Number of eruptions');

plt.plot(x_plot,sp.stats.beta.pdf(x_plot,0.4694039854411987, 1.7452248296341))
plt.plot(x_plot,sp.stats.beta.pdf(x_plot,1, 1))

post_pdf_contribs = sp.stats.beta.pdf(np.atleast_3d(x_plot),
                                      trace['mean'][:, np.newaxis, :],
                                      [:, np.newaxis, :])
post_pdfs = (trace['w'][:, np.newaxis, :] * post_pdf_contribs).sum(axis=-1)

post_pdf_low, post_pdf_high = np.percentile(post_pdfs, [2.5, 97.5], axis=0)

help(pm.NUTS)

help(pm.MvNormal)

with model:
    start = pm.find_MAP()
    print(start)
start


with model:
    start = pm.find_MAP()

    step1 = pm.NUTS([p, sigma, mean])
    step2 = pm.BinaryMetropolis([component], scaling=.01)
    trace = pm.sample(240000, start=start, step=[step1, step2])
    pm.traceplot(trace[40000:][::50], varnames=["p","sigma","mean"])
    print([int(model.logp(t)) for t in trace[40000:][::400]])
    print([trace[i]['component'].sum() for i in np.arange(0,len(trace),400)])

pm.traceplot(trace[200000:][::50], varnames=["p","sigma","mean"])



pm.traceplot(trace[:][::50], varnames=["p","sigma","mean"])

sns.distplot(data)


with pm.Model() as model:
    sigma = pm.Lognormal('sigma', np.zeros(2), np.ones(2), shape=2)

    nu = pm.Uniform('nu', 0, 5)
    C_triu = pm.LKJCorr('C_triu', nu, 2)

with model:
    C = pm.Deterministic('C', T.fill_diagonal(C_triu[np.zeros((2, 2), dtype=np.int64)], 1.))

    sigma_diag = pm.Deterministic('sigma_mat', T.nlinalg.diag(sigma))
    cov = pm.Deterministic('cov', T.nlinalg.matrix_dot(sigma_diag, C, sigma_diag))
    tau = pm.Deterministic('tau', T.nlinalg.matrix_inverse(cov))

with model:
    mu = pm.Uniform('mean', 0., 1., shape=K)

    x_ = pm.MvNormal('x', mu, tau, observed=x)
