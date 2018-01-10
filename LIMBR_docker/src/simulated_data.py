%matplotlib inline
from random import *
import numpy as np
import seaborn as sns
import pandas as pd
import os
from numpy.linalg import svd
from fbpca import pca
import time
from sklearn.preprocessing import scale, MinMaxScaler
from sklearn.decomposition import RandomizedPCA
from scipy.optimize import curve_fit
from pylab import *
from math import atan2
from scipy.stats import linregress
import itertools
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.metrics import roc_curve, auc

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


exp_class = [j for i in range(0,12) for j in [i]*3]*2

classes = exp_class

len(classes)



classes

circ = [randint(0,1) for b in range(1,10001)]

circ

base = np.arange(0,(4*np.pi),(4*np.pi/24))


sim = []
phases = []

for i in circ:
    if i == 1:
        temp=[]
        p = randint(0,1)
        phases.append(p)
        temp.append(np.sin(base+np.random.normal(0,0.4,1)+np.pi*p)+np.random.normal(0,1.75,24))
        temp.append(np.sin(base+np.random.normal(0,0.4,1)+np.pi*p)+np.random.normal(0,1.75,24))
        temp.append(np.sin(base+np.random.normal(0,0.4,1)+np.pi*p)+np.random.normal(0,1.75,24))
        temp2 = []
        for i in range(len(temp[0])):
            temp2.append(temp[0][i])
            temp2.append(temp[1][i])
            temp2.append(temp[2][i])
        sim.append(temp2)
    else:
        phases.append('nan')
        sim.append(np.random.normal(0,1,72))

np.shape(sim)
simdf = pd.DataFrame(sim,columns=csp_norm.columns.values)
simdf.index.names = ['#']

simnoise = []
trend1list = []
#trend2list = []
trend3list = []
for i in sim:
    temp = []
    t1 = randint(0,1)
    #t2 = randint(0,1)
    t3 = randint(0,1)
    trend1list.append(t1)
    #trend2list.append(t2)
    trend3list.append(t3)
    for j in range(len(i)):
        trend = [i*t1 for i in ([0,0,3]*24)]
        #trend2 = [k*t2*.5 for k in ([0]*63 + [3]*9)]
        trend3 = [i*t3 for i in ([0,3,0]*12) + [0,0,0]*12]
        temp.append(i[j]+trend[j]+trend3[j])
        #temp.append(i[j]+trend[j])
    simnoise.append(temp)


simndf = pd.DataFrame(simnoise,columns=csp_norm.columns.values)
simndf.index.names = ['#']


simdf.to_csv('simulated_data_baseline.txt',sep='\t')
simndf.to_csv('simulated_data_with_noise.txt',sep='\t')
k = pd.concat([pd.Series(circ),pd.Series(phases),pd.Series(trend1list),pd.Series(trend3list)],axis=1)
#k = pd.concat([pd.Series(circ),pd.Series(phases),pd.Series(trend1list)],axis=1)
k.columns = ['circ','phase','trend1','trend3']
#k.columns = ['circ','phase','trend1']

k.to_csv('simulated_data_key.txt',sep='\t')


#simdf= pd.read_csv('simulated_data_baseline.txt',sep='\t')
#simdf = simdf.set_index('#')

sns.clustermap(simdf.head(n=100),col_cluster=False)
plt.savefig('simulated_data.pdf')



sns.clustermap(simndf.head(n=100),col_cluster=False)
plt.savefig('simulated_noised_data.pdf')

len(tpoints)

cors = prim_cor(simndf.values,tpoints,12)

cors = pd.Series(cors,index=simndf.index)

merged = pd.concat([cors, k['circ'].astype('bool')], axis=1,join='inner')
merged.columns = ['cors','circ']

sns.violinplot(x="circ", y="cors", data=merged)
plt.title('Correlation in Simulated Data with Noise')
plt.savefig('primary_corellations_simulated_noised.pdf')

cors = prim_cor(simdf.values,tpoints,12)

cors = pd.Series(cors,index=simdf.index)

merged = pd.concat([cors, k['circ'].astype('bool')], axis=1,join='inner')
merged.columns = ['cors','circ']

sns.violinplot(x="circ", y="cors", data=merged)
plt.title('Correlation in Simulated Data without Noise')
plt.savefig('primary_corellations_simulated.pdf')

simj = pd.read_csv('simulated_data_baseline__jtkout_GammaP.txt',sep='\t')
simj = simj.set_index('ID')

sns.distplot(simj['GammaP'])

merged = []
merged = pd.concat([simj['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)
plt.title('Initial Data')
plt.savefig('simulated_classification.pdf')

merged['circ'].values
-merged['p'].values

fpr, tpr, _ = roc_curve(merged['circ'].values, -merged['p'].values)
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, label='ROC curve initial data (area = %0.2f)' % roc_auc)
plt.plot(fprdn, tprdn, label='ROC curve denoised Lowess (area = %0.2f)' % roc_aucdn)
plt.plot(fprdncr, tprdncr, label='ROC curve denoised circadian replicate (area = %0.2f)' % roc_aucdncr)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('ROC_curves.pdf')

simnj = pd.read_csv('simulated_data_with_noise__jtkout_GammaP.txt',sep='\t')
simnj = simnj.set_index('ID')


sns.distplot(simnj['GammaP'])


merged = pd.concat([simnj['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)
plt.title('Data with Noise')
plt.savefig('simulated_classification_noised.pdf')

fprn, tprn, _ = roc_curve(merged['circ'].astype('bool').values, -merged['p'].values)
roc_aucn = auc(fprn, tprn)


cors = prim_cor(simndf.values,tpoints,12)

uncor = [(i<(np.percentile(cors,25))) for i in cors]

simndf_reduced = simndf[uncor]

def get_res(arr,l):
    res = []
    for row in arr:
        ys = lowess(row, l,delta=4)[:,1]
        res.append(row - ys)
    return res

res = get_res(simndf_reduced.values,tpoints)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,tpoints,100,tks_)
end = time.time()
print(end - start)

sigs

ps = eig_reg(simndf_reduced.values,res,sigs,.05)

ts = subset_svd(simndf_reduced.values,ps,0.5)

fin_res = np.dot(linalg.lstsq(np.asarray(ts).T,simndf.values.T)[0].T,np.asarray(ts))

svd_norm = simndf.values - fin_res

svd_norm = pd.DataFrame(svd_norm,index=simndf.index,columns=simndf.columns)
svd_norm = pd.DataFrame(scale(svd_norm.values,axis=1),columns=svd_norm.columns,index=svd_norm.index)

svd_norm.index.names = ['#']

svd_norm.to_csv('simulated_data_denoised.txt',sep='\t')

svd_norm = pd.read_csv('simulated_data_denoised.txt',sep='\t')
svd_norm = svd_norm.set_index('#')

sns.clustermap(svd_norm.head(n=100),col_cluster=False)
plt.savefig('simulated_denoised_data.pdf')

simdnj = pd.read_csv('simulated_data_denoised__jtkout_GammaP.txt',sep='\t')
simdnj = simdnj.set_index('ID')


sns.distplot(simdnj['GammaP'])


merged = pd.concat([simdnj['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)
plt.title('Lowess Denoised')
plt.savefig('simulated_classification_denoised.pdf')

fprdn, tprdn, _ = roc_curve(merged['circ'].astype('bool').values, -merged['p'].values)
roc_aucdn = auc(fprdn, tprdn)




svd_norm = normalize(simndf.values,res,sigs,0.05)


svd_norm = pd.DataFrame(svd_norm,index=simndf.index,columns=simndf.columns)

svd_norm.index.names = ['#']

svd_norm.to_csv('simulated_denoised_circ_rep.txt',sep='\t')

sns.clustermap(svd_norm.head(n=100),col_cluster=False)

simdnj = pd.read_csv('simulated_denoised_circ_rep__jtkout_GammaP.txt',sep='\t')
simdnj = simdnj.set_index('ID')


sns.distplot(simdnj['GammaP'])

merged = pd.concat([simdnj['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)

exp_class = [j for i in range(0,24) for j in [i]*3]

classes = exp_class

len(classes)

classes


def get_res(arr,l):
    m = {}
    for v in set(l):
        indices = [i for i, x in enumerate(l) if x == v]
        m[v] = np.mean(arr[:,indices],axis=1)
    ma = np.zeros(np.shape(arr))
    for i in range(len(l)):
        ma[:,i]=m[l[i]]
    return np.subtract(arr,ma)

res = get_res(simndf.values,classes)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,classes,100,tks_)
end = time.time()
print(end - start)

sigs

svd_norm = normalize(simndf.values,res,sigs,0.05)


svd_norm = pd.DataFrame(svd_norm,index=simndf.index,columns=simndf.columns)

svd_norm.index.names = ['#']

svd_norm.to_csv('simulated_data_denoised_rep.txt',sep='\t')

simdnjcr = pd.read_csv('simulated_data_denoised_circ_rep__jtkout_GammaP.txt',sep='\t')
simdnjcr = simdnjcr.set_index('ID')


sns.distplot(simdnjcr['GammaP'])


merged = pd.concat([simdnjcr['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)
plt.title('Circadian Replicate Denoised')
plt.savefig('simulated_classification_denoised_circrep.pdf')

fprdncr, tprdncr, _ = roc_curve(merged['circ'].astype('bool').values, -merged['p'].values)
roc_aucdncr = auc(fprdncr, tprdncr)

simdnjr = pd.read_csv('simulated_data_denoised_rep__jtkout_GammaP.txt',sep='\t')
simdnjr = simdnjr.set_index('ID')


sns.distplot(simdnjr['GammaP'])


merged = pd.concat([simdnjr['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)

fprdnr, tprdnr, _ = roc_curve(~merged['circ'].astype('bool').values, merged['p'].values)
roc_aucdnr = auc(fprdnr, tprdnr)

sns.clustermap(svd_norm.head(n=100),col_cluster=False)

simdnj2 = pd.read_csv('simulated_denoised_rep__jtkout_GammaP.txt',sep='\t')
simdnj2 = simdnj2.set_index('ID')


sns.distplot(simdnj2['GammaP'])

merged = pd.concat([simdnj2['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)

merged = pd.concat([simdnjcr['Phase'], k['circ'],k['phase']], axis=1,join='inner')
merged.columns = ['Phase_pred','circ','Phase_real']

sns.jointplot(x='Phase_pred',y='Phase_real',data=merged[merged['circ']==1],kind='kde')

merged = pd.concat([simdnj['Phase'], k['circ'],k['phase']], axis=1,join='inner')
merged.columns = ['Phase_pred','circ','Phase_real']

phasediff = simj['Phase'] - simdnj['Phase']

theta = np.linspace(0,2*np.pi,145)
theta_matrix = np.meshgrid(theta)[0].tolist()
transformed = [(i % 23)*np.pi*2/23 for i in phasediff.values.tolist()]

transformed = transformed + [(i - (2*np.pi)) for i in transformed] + [(i + (2*np.pi)) for i in transformed]

kernel = stats.gaussian_kde(transformed + [(i + (2*np.pi)) for i in transformed])
Z = kernel(theta_matrix)

kde = KernelDensity(kernel='gaussian',bandwidth=0.25).fit(np.asarray(transformed).reshape(-1,1))
Z = kde.score_samples(np.asarray(theta_matrix).reshape(-1,1))

Z  = scale(Z)
Z = Z - np.min(Z)
min_max_scaler = MinMaxScaler()
Z = min_max_scaler.fit_transform(Z.reshape(-1,1))

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(70)
CS = plt.plot(np.linspace(0,2*np.pi,145), Z, c= 'b')
plt.savefig('simulated_phasediffs.pdf')


sns.jointplot(x='Phase_pred',y='Phase_real',data=merged[merged['circ']==1],kind='kde')

merged = pd.concat([simj['Phase'], k['circ'],k['phase']], axis=1,join='inner')
merged.columns = ['Phase_pred','circ','Phase_real']

sns.jointplot(x='Phase_pred',y='Phase_real',data=merged[merged['circ']==1],kind='kde')
plt.savefig('initial_phases.pdf')

def get_res(arr,l):
    res = []
    for row in arr:
        def harm(x, p1,p2,p3):
            return p1*np.cos(2*np.pi*p2*x + 2*np.pi*p3)
        amplitude = row.max() - row.min()
        popt, pcov = curve_fit(harm, l, row, p0=(amplitude,.043478261,0), bounds=([0,.043478260,-np.inf], [np.inf, .043478262, np.inf]))
        res.append(row - harm(l,popt[0],popt[1],popt[2]))
    return res

def get_tpoints(l):
    tpoints = [i.replace('CT','') for i in l]
    tpoints = [int(i.split('_')[0]) for i in tpoints]
    return np.asarray(tpoints)


tpoints = get_tpoints(simndf.columns.values)

res = get_res(simndf.values,tpoints)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,tpoints,100,tks_)
end = time.time()
print(end - start)

sigs


def autocorr(l,shift):
    return dot(l, np.roll(l, shift)) / dot(l, l)


def prim_cor(arr,l,per,reps):
    cors = []
    for row in arr:
        ave = [np.mean(i) for i in [row[j:j + reps] for j in range(0, len(row), reps)]]
        cors.append((autocorr(ave,per) - autocorr(ave,(per//2))))
    return cors

temp = [np.mean(i) for i in [simdf.values[1][j:j + 3] for j in range(0, len(simdf.values[1]), 3)]]

len(temp)
autocorr(temp,12)
autocorr(temp,6)

cors = prim_cor(simndf.values,tpoints,12)

plt.scatter(k['circ'].values.tolist(),cors)

sns.violinplot(x=k['circ'].values.tolist(),y=cors)

k['circ'].values.tolist()

np.mean(cors)
uncor = [(i<(np.mean(cors)*.25)) for i in cors]
len(simndf[uncor])
simndf_reduced = simndf[uncor]

tpoints = get_tpoints(simndf_reduced.columns.values)
#res = get_res(simndf_reduced.values,classes)
res = get_res(simndf_reduced.values,tpoints)

sns.clustermap(res[:1000],col_cluster=False)

tks_ = get_tks(res)
np.sum(tks_)

start = time.time()
sigs = perm_test(res,tpoints,100,tks_)
end = time.time()
print(end - start)

sigs

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

def est_pi_naught(probs,lam):
    return len([i for i in probs if i > lam])/(len(probs)*lam)

def est_pi_sig(probs,l):
    pi_0 = est_pi_naught(probs,l)
    sp = np.sort(probs)
    return sp[int(floor(pi_0*len(probs)))]

def subset_svd(arr,plist,lam):
    trends = []
    for entry in plist:
        sub = []
        thresh = est_pi_sig(entry,lam)
        for i in range(len(entry)):
            if entry[i] < thresh:
                sub.append(arr[i])
        U, s, V = np.linalg.svd(sub)
        trends.append(V.T[:,0])
    return trends


ps = eig_reg(simndf_reduced.values,res,sigs,.05)
sns.distplot(ps[0])
est_pi_naught(ps[0],.5)

est_pi_sig(ps[0],.5)

ts = subset_svd(simndf_reduced.values,ps,0.5)

plt.plot(range(len(ts[0])),ts[0])

plt.plot(range(len(ts[1])),ts[1])

plt.plot(range(len(ts[2])),ts[2])

U, s, V = np.linalg.svd(res)
sig = V.T[:,:len([i for i in itertools.takewhile(lambda x: x < .05, sigs)])]


np.shape(ts)
np.shape(sig)
plt.plot(range(len(sig.T[0])),sig.T[0])

plt.plot(range(len(sig.T[1])),sig.T[1])

plt.plot(range(len(sig.T[2])),sig.T[2])

np.shape(simndf_reduced.values)

#fin_res = np.dot(linalg.lstsq(np.asarray(sig),simndf.values.T)[0].T,np.asarray(sig).T)
fin_res = np.dot(linalg.lstsq(np.asarray(ts).T,simndf.values.T)[0].T,np.asarray(ts))

sns.clustermap(fin_res[:100],col_cluster=False)

svd_norm = simndf.values - fin_res


svd_norm = pd.DataFrame(svd_norm,index=simndf.index,columns=simndf.columns)

svd_norm.index.names = ['#']

svd_norm.to_csv('simulated_denoised_new_lowess.txt',sep='\t')

sns.clustermap(svd_norm.head(n=1000),col_cluster=False)

simdnj = pd.read_csv('simulated_denoised_new__jtkout_GammaP.txt',sep='\t')
simdnj = simdnj.set_index('ID')


sns.distplot(simdnj['GammaP'])

merged = pd.concat([simdnj['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)

simdnj = pd.read_csv('simulated_denoised_new2__jtkout_GammaP.txt',sep='\t')
simdnj = simdnj.set_index('ID')


sns.distplot(simdnj['GammaP'])

merged = pd.concat([simdnj['GammaP'], k['circ']], axis=1,join='inner')
merged.columns = ['p','circ']

merged[merged['circ']==1]['p'].mean()

len(merged[(merged['circ']==1) & (merged['p']<.05)])
len(merged[(merged['circ']==0) & (merged['p']<.05)])

sns.violinplot(x="circ", y="p", data=merged)

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
