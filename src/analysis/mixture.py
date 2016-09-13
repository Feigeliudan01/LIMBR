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
theano.config.compute_test_value = 'warn'
theano.config.compute_test_value = 'off'
import theano
import time
import os

os.chdir('/Users/user4574/Dropbox (Dunlap-Loros labs)/MSLearn')

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
    trace = pm.sample(200, step=step, start= start)


with model:
    start = pm.find_MAP()
    trace = pm.sample(200, start=start)

pm.traceplot(trace, varnames=["a"])



#takes ~18min with 19x19 grid
np.mean(trace['a'][::50,0])
N= 1000
delta = 0.05
xgrid = np.arange(delta, 1.0, delta)
ygrid = np.arange(delta, 1.0, delta)
X, Y = np.meshgrid(xgrid, ygrid)
np.shape(X)
a = [84.193916952007285,0.16078286199088215,0.20358242471646312,0.1276771579030774]
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
plt.savefig('woooo.pdf')

#use switch and draw from either bivariate beta, beta on one axis random uniform on the other or random uniform on both

sns.jointplot(x=np.reshape(bwprot['GammaP'].values[:1000],(1,-1)), y = np.reshape(brna['GammaP'].values[:1000],(1,-1)), kind='kde')
sns.jointplot(x=x, y = y)
##############################################################################################################################
#mixture

ndata = len(x[0])
with pm.Model() as model:
    a = pm.Gamma('a',10.,10.,shape=[4,4])
    #a1 = pm.Gamma('a1',1,1)
    #a2 = pm.Gamma('a2',1,1)
    #b1 = pm.Gamma('b1',1,1)
    #b2 = pm.Gamma('b2',1,1)
    pun = pm.Uniform('pun',0,1)
    dd = pm.Dirichlet('dd', a=np.array([1., 1., 1.,1.]), shape=4)
    category = pm.Categorical(name='category', p=prob_dist, shape=ndata)
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
    def xonly(x,y):
        return scipy.stats.beta(a1, b1).pdf(x)
    def yonly(x,y):
        return scipy.stats.beta(a2, b2).pdf(y)
    def un(x,y):
        return pun
    funcs = np.asarray([dens,xonly])
    funcs = theano.shared(funcs)
    #real_cat = T.switch(2 > category, 0, 1)
    #funclist = T.vector()
    #funcfunc = theano.function([real_cat,funclist], funclist[real_cat])
    Bibeta = pm.DensityDist("Bibeta", dens,  observed={"x": x, 'y':y})

ndata = len(x[0])
with pm.Model() as model:
    #alist = pm.Gamma('a',mu=1,sd=1,shape=(4,4))
    a1 = pm.Gamma('a1',1.,1.,shape=2)
    a2 = pm.Gamma('a2',1.,1.,shape=2)
    a3 = pm.Gamma('a3',1.,1.,shape=2)
    a4 = pm.Gamma('a4',1.,1.,shape=2)
    dd = pm.Dirichlet('dd', a=np.array([1., 1.]), shape=2)
    category = pm.Categorical(name='category', p=dd,shape=ndata)
    #a = alist[category]
    a1_r = a1[category]
    a2_r = a2[category]
    a3_r = a3[category]
    a4_r = a4[category]
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
        #need to fix next line so that both u and as are getting iterated over seperately
        myint = lambda u: ((u**(a4_r-1))*((x-u)**(a3_r-1))*((y-u)**(a2_r-1))*((1-x-y+u)**(a1_r-1)))
        num = T.dscalar('num')
        ma = T.dscalar('ma')
        mi = T.dscalar('mi')
        tlinspace = lambda ma,mi,num: T.dot(((T.reshape(T.arange(num),(-1,1)) + .5) / num),(ma-mi))+T.extra_ops.repeat(mi,num,axis=0)
        lin = tlinspace(up,low,N)
        area = T.sum(myint(lin))*(up-low)/N
        return T.log(T.gamma(T.sum([a1_r,a2_r,a3_r,a4_r]))/T.prod(T.gamma([a1_r,a2_r,a3_r,a4_r]))*area)
    Bibeta = pm.DensityDist("Bibeta", dens,  observed={"x": x, 'y':y})

model

# setup model
model = pm.Model()
k = 3
with model:
    # cluster sizes
    a = pm.constant(np.array([1., 1., 1.]))
    p = pm.Dirichlet('p', a=a, shape=k)
    # ensure all clusters have some points
    p_min_potential = pm.Potential('p_min_potential', tt.switch(tt.min(p) < .1, -np.inf, 0))


    # cluster centers
    means = pm.Normal('means', mu=[0, 0, 0], sd=15, shape=k)
    # break symmetry
    order_means_potential = pm.Potential('order_means_potential',
                                         tt.switch(means[1]-means[0] < 0, -np.inf, 0)
                                         + tt.switch(means[2]-means[1] < 0, -np.inf, 0))

    # measurement error
    sd = pm.Uniform('sd', lower=0, upper=20)

    # latent cluster of each observation
    category = pm.Categorical('category',
                           p=p,
                           shape=ndata)

    # likelihood for each observed value
    points = pm.Normal('obs',
                     mu=means[category],
                     sd=sd,
                     observed=x)

# fit model
with model:
    step1 = pm.Metropolis(vars=[p])
    step2 = pm.Metropolis(vars=[sd, means])
    step3 = SequentialScanDiscreteMetropolis(var=category, values=[0,1,2])
    tr = pm.sample(10000, step=[step1, step2] + [step3]*ndata)

with model:
    start = pm.find_MAP()
    start

start

with model:
    start = pm.find_MAP()
    start['a1_log'] = np.array([1.,1.])
    start['a2_log'] = np.array([1.,1.])
    start['a3_log'] = np.array([1.,1.])
    start['a4_log'] = np.array([1.,1.])
    trace = pm.sample(200, start=start)

pm.plots.traceplot(trace)

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

with model:
    step1 = pm.Metropolis(vars=[dd])
    step2 = pm.Metropolis(vars=[a1,a2,a3,a4])
    step3 = SequentialScanDiscreteMetropolis(var=category, values=[0,1])
    tr = pm.sample(200, step=[step1, step2, step3])

with model:
    #start = pm.find_MAP()
    trace = pm.sample(10)

##############################################################################################################################
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
