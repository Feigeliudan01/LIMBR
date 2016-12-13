import pymc3 as pm, theano.tensor as tt
from pymc3.math import switch
import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from pymc3.backends import SQLite
from scipy.stats import beta
import sys
import optparse
from tqdm import tqdm
from scipy.integrate import quad
from scipy.optimize import brentq
import pickle

def threshold(pvals,fdr):
    k=2
    y = pvals
    model = pm.Model()
    with pm.Model() as model:
        # cluster sizes
        p = pm.Dirichlet('p', a=np.array([1., 1.]), shape=k)
        # ensure all clusters have some points
        p_min_potential = pm.Potential('p_min_potential', switch(tt.min(p) < .1, -np.inf, 0))
        # latent cluster of each observation
        category = pm.Categorical('category',p=p,shape=len(y))
        a = pm.HalfCauchy('a', beta=[1,1],shape=k)
        b = pm.HalfCauchy('b', beta=[1,1],shape=k)
        points = pm.Beta('obs',alpha=a[category], beta=b[category], observed=y)
    # fit model
    with model:
        step1 = pm.Metropolis(vars=[p, a, b])
        step2 = pm.ElemwiseCategorical(vars=[category], values=[0, 1])
        tr = pm.sample(10000, step=[step1, step2])
    w = np.mean(tr['category',1000::5])
    aest1 = np.mean(tr['a',1000::5][:,0])
    best1 = np.mean(tr['b',1000::5][:,0])
    aest2 = np.mean(tr['a',1000::5][:,1])
    best2 = np.mean(tr['b',1000::5][:,1])
    return [w, aest1, best1, aest2, best2]

    if aest1/(aest1+best1) > aest2/(aest2+best2):
        rv2 = beta(aest1, best1)
        rv1 = beta(aest2,best2)
        def integrand(p):
            return (rv2.pdf(p)*(1-w))/((w*rv1.pdf(p))+((1-w)*rv2.pdf(p)))
    else:
        rv1 = beta(aest1, best1)
        rv2 = beta(aest2,best2)
        def integrand(p):
            return (rv2.pdf(p)*w)/(((1-w)*rv1.pdf(p))+(w*rv2.pdf(p)))
    def func(x):
        y, err = quad(integrand, 0, x)
        return y - fdr
    return brentq(func, 0, 1)


thresh = {}
n = 0
for i in list(set(net_data.index.get_level_values('experiment'))):
    n += 1
    print(n)
    pvals = net_data[net_data.index.get_level_values('experiment') == i]['p_value'].values.tolist()
    thresh[i] = threshold(pvals,10**(-9))
    pickle.dump( thresh, open( "thresholds.p", "wb" ) )


def gen_threshold(l, fdr):
        w, aest1, best1, aest2, best2 = l
        if aest1/(aest1+best1) > aest2/(aest2+best2):
            rv2 = beta(aest1, best1)
            rv1 = beta(aest2,best2)
            def integrand(p):
                return (rv2.pdf(p)*(1-w))/((w*rv1.pdf(p))+((1-w)*rv2.pdf(p)))
        else:
            rv1 = beta(aest1, best1)
            rv2 = beta(aest2,best2)
            def integrand(p):
                return (rv2.pdf(p)*w)/(((1-w)*rv1.pdf(p))+(w*rv2.pdf(p)))
        def func(x):
            y, err = quad(integrand, 0, x)
            return y - fdr
        return brentq(func, 0, 1)

def gen_thresh_mat(d,f):
    newd = {}
    for key, value in d.iteritems():
        newd[key] = gen_threshold(value,f)
    return newd

def main(argv):
    parser = optparse.OptionParser()
    parser.add_option('-i','--input', dest='inputfile')
    parser.add_option('-p','--precision', dest='precision')
    (options, args) = parser.parse_args()
    data = pd.read_csv(options.inputfile,sep='\t',header=None)
    y = data[0].values

def gen_network(df,tdict):
    newdf = pd.DataFrame(np.zeros((len(df),1)),index=df.index,columns=['binding'])
    for i in list(set(net_data.index.get_level_values('experiment'))):
        newdf.ix[newdf.index.get_level_values('experiment') == i,'binding'] = net_data.ix[net_data.index.get_level_values('experiment') == i,'p_value'] < tdict[i]
    return newdf

def gen_elbow(df,bs,minfdr_exp,maxfdr_exp):
    l = []
    for i in range(minfdr_exp,maxfdr_exp):
        tmat = gen_thresh_mat(bs,(10**-i))
        net = gen_network(df,tmat)
        l.append(net['binding'].groupby(net.index.get_level_values('chip_protein')).sum().values.tolist())
        #l.append(net.groupby(net.index.get_level_values('target_gene')).max().sum().values.tolist()[0])
        #l.append(net['binding'].sum())
    names = net.index.values.tolist()
    return l, names

def gen_elbow(df,bs,minfdr_exp,maxfdr_exp):
    l = []
    for i in np.linspace(minfdr_exp,maxfdr_exp,10):
        tmat = gen_thresh_mat(bs,i)
        net = gen_network(df,tmat)
        #l.append(net['binding'].groupby(net.index.get_level_values('chip_protein')).sum().values.tolist())
        l.append(net.groupby(net.index.get_level_values('target_gene')).max().sum().values.tolist()[0])
        #l.append(net['binding'].sum())
    names = net.index.values.tolist()
    return l, names

pab = pickle.load( open( "thresholds.p", "rb" ) )

sig = gen_network(net_data,gen_thresh_mat(pab,10**-6))

#sig = sig[sig['binding'] == True]

conditions = pd.read_csv('TF_exp_codes.txt',sep='\t')

sig_con = pd.merge(sig.reset_index(), conditions[['real_con','File.Name']], left_on='experiment', right_on='File.Name')

sig_con = sig_con.drop('File.Name',axis=1)

def gen_bmat(df):
    tfs = list(set(df['chip_protein']))
    temp = pd.DataFrame(np.zeros((len(tfs),len(tfs))),index=tfs,columns=tfs)
    for i in tfs:
        pos = df[(df['chip_protein']==i) & (df['binding']==True)][['target_gene','real_con']]
        pos = [list(x) for x in pos.values]
        for j in tfs:
            sub = df[df['chip_protein'] == j]
            indeces = []
            for k in pos:
                indeces.append(pair_lookup(sub,k[0],k[1]))
            match = sub[np.asarray(indeces).any(axis=0)]
            if len(match)>0:
                temp = temp.set_value(i,j, float(len(match[match['binding']==True])/len(match)))
            else:
                temp = temp.set_value(i,j, 0)
    return temp


        for j in list(set(df['real_con'])):

        for j in list(set(df['target_gene'])):



def pair_lookup(df,tgene,con):
    return (df['target_gene']==tgene) & (df['real_con']==con)

targets = set(sig.index.get_level_values('target_gene').values.tolist())
tfs = set(sig.index.get_level_values('chip_protein').values.tolist())

bmat = pd.DataFrame(np.zeros((len(targets),len(tfs))),index=targets,columns=tfs)

pairs = sig.reset_index().ix[:,:2]

for i in range(len(pairs)):
    bmat = bmat.set_value(pairs.ix[i,1],pairs.ix[i,0], 1)

bmat.to_csv('binary_occupancy_matrix.txt',sep='\t')

net_data[]

    newpvals = []
    for i in tqdm(data[0].values):
            newpvals.append(FDR(i,int(options.precision)))
            #newpvals.append(FDR2(i))
            #newpvals.append(FDR3(i))

    final = pd.DataFrame(np.zeros((len(data.index),len(data.columns))),index=data.index,columns=data.columns)
    final[0] = newpvals
    final[1] = data[1]
    final.to_csv(options.inputfile[:-4]+'_adjusted.txt',sep='\t',header=False,index=False)



if __name__ == '__main__':
    main(sys.argv[1:])


#pm.plots.traceplot(tr, ['p', 'sd', 'means']);
#pm.plots.traceplot(tr[5000::5], ['p', 'a1', 'b1']);
#x = np.linspace(0.001,0.999, 1000)
#plt.plot(x, beta.pdf(x, aest, best),'r-', alpha=0.6, label='beta pdf')
#plt.show()
