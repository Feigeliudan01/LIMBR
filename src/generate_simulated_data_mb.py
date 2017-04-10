from random import *
import numpy as np
import pandas as pd
import pickle
from sklearn.preprocessing import scale

seed(4574)
cols = pickle.load( open( 'src/labels.p', "rb" ) )
cols = cols[2:]
randBinList = lambda n: [randint(0,1) for b in range(n)]
circ = [randint(0,1) for b in range(10000)]

base = np.arange(0,(4*np.pi),(4*np.pi/24))

def gen_sim_data(suf):
    sim = []
    phases = []

    for i in circ:
        if i == 1:
            temp=[]
            p = randint(0,1)
            phases.append(p)
            temp.append(np.sin(base+np.random.normal(0,0.25,1)+np.pi*p)+np.random.normal(0,1,24))
            temp.append(np.sin(base+np.random.normal(0,0.25,1)+np.pi*p)+np.random.normal(0,1,24))
            temp.append(np.sin(base+np.random.normal(0,0.25,1)+np.pi*p)+np.random.normal(0,1,24))
            temp2 = []
            for i in range(len(temp[0])):
                temp2.append(temp[0][i])
                temp2.append(temp[1][i])
                temp2.append(temp[2][i])
            sim.append(temp2)
        else:
            phases.append('nan')
            sim.append(np.random.normal(0,1,72))

    simdf = pd.DataFrame(np.asarray(sim),columns=cols)
    simdf.index.names = ['#']
    simdf = pd.DataFrame(scale(simdf.values,axis=1),columns=simdf.columns,index=simdf.index)

    simnoise = []
    trend1list = []
    trend2list = []
    trend3list = []
    basetrend = np.random.normal(0.,2.,72)
    basetrend2 = np.random.normal(0.,2.,72)
    basetrend3 = np.random.normal(0.,2.,72)
    for i in sim:
        temp = []
        t1 = randint(0,1)
        t2 = randint(0,1)
        t3 = randint(0,1)
        trend1list.append(t1)
        trend2list.append(t2)
        trend3list.append(t3)
        for j in range(len(i)):
            trend = [k*t1 for k in basetrend]
            trend2 = [k*t2 for k in basetrend2]
            trend3 = [k*t3 for k in basetrend3]
            temp.append(i[j]+trend[j]+trend2[j]+trend3[j])
        simnoise.append(temp)

    simndf = pd.DataFrame(np.asarray(simnoise),columns=cols)
    simndf.index.names = ['#']
    simndf = pd.DataFrame(scale(simndf.values,axis=1),columns=simndf.columns,index=simndf.index)

    simdf.to_csv('output/simdata/simulated_data_baseline_'+str(m)+'.txt',sep='\t')
    simndf.to_csv('output/simdata/simulated_data_with_noise_'+str(m)+'.txt',sep='\t')
    simndf.insert(0, 'Peptide', ['1']*len(simndf))
    simndf.index.names = ['Protein']
    simndf.to_csv('output/simdata/simulated_data_with_noise_for_sva_'+str(m)+'.txt',sep='\t')
    #k = pd.concat([pd.Series(circ),pd.Series(phases),pd.Series(trend1list),pd.Series(trend2list)],axis=1)
    k = pd.concat([pd.Series(circ),pd.Series(phases),pd.Series(trend1list),pd.Series(trend2list),pd.Series(trend3list)],axis=1)
    k.columns = ['circ','phase','trend1','trend2','trend3']
    k.to_csv('output/simdata/simulated_data_key_'+str(m)+'.txt',sep='\t')
    trends=pd.DataFrame([basetrend,basetrend2,basetrend3])
    trends.to_csv('output/simdata/trends_'+str(m)+'.txt',sep='\t')

for m in range(1,6):
    gen_sim_data(m)
