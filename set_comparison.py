%matplotlib inline
import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats
from sklearn.preprocessing import normalize, MinMaxScaler
import matplotlib.cm as cm
from scipy.stats import beta

nwprot = pd.read_csv('wt_lowess__jtkout_GammaP.txt',sep='\t')
ncprot = pd.read_csv('csp_lowess__jtkout_GammaP.txt',sep='\t')
repwprot = pd.read_csv('rep_wt_svd_norm__jtkout_GammaP.txt',sep='\t')
reprna =  pd.read_csv('rep_rna_svd_norm__jtkout_GammaP.txt',sep='\t')
oldwprot = pd.read_csv('JTK.wt_eigenms_imputed_old.txt',sep='\t')
midwprot = pd.read_csv('wt_SVD_normalized__jtkout_GammaP.txt',sep='\t')
circrepwprot = pd.read_csv('circ_rep_wt_svd_norm__jtkout_GammaP.txt',sep='\t')
circrepsepwprot = pd.read_csv('circ_rep_sep_wt_svd_norm__jtkout_GammaP.txt',sep='\t')

wprot_in = pd.read_csv('wt_reduced_subset_lowess_qnorm_final.txt',sep='\t')
wprot_in = wprot_in.set_index('#')

wprot_in = pd.DataFrame(scale(wprot_in.values,axis=1),columns=wprot_in.columns,index=wprot_in.index)

bwprot[bwprot['GammaP']<.05][0]

len(bwprot)
len(bcprot)

wprot_in.ix[bwprot[bwprot['GammaP']<.05].index]

sns.clustermap(wprot_in.ix[bwprot[bwprot['GammaP']<.05].index].head(n=1000),col_cluster=False)
plt.savefig('wt_circ_1000.pdf')

cprot_res = pd.read_csv('csp_residuals__jtkout_GammaP.txt',sep='\t')
cprot_res = cprot_res.set_index('ID')

wprot_res = pd.read_csv('wt_residuals__jtkout_GammaP.txt',sep='\t')
wprot_res = wprot_res.set_index('ID')

brna = pd.read_csv('rna_reduced_subset_lowess_qnorm_final_cor__jtkout_GammaP.txt',sep='\t')
brna = brna.set_index('ID')

bcprot = pd.read_csv('csp_reduced_subset_lowess_qnorm_final_cor__jtkout_GammaP.txt',sep='\t')
bcprot = bcprot.set_index('ID')


bwprot = pd.read_csv('wt_reduced_subset_lowess_qnorm_final_cor__jtkout_GammaP.txt',sep='\t')
bwprot = bwprot.set_index('ID')

bwprot.ix['NCU09482']

nwprot = nwprot.set_index('ID')
ncprot = ncprot.set_index('ID')
repwprot = repwprot.set_index('ID')
reprna = reprna.set_index('ID')
midwprot = midwprot.set_index('ID')
circrepwprot = circrepwprot.set_index('ID')
circrepsepwprot = circrepsepwprot.set_index('ID')

oldwprot = oldwprot.set_index('rownames(data)')
oldwprot.index.names = ['ID']
oldwprot = oldwprot.groupby(oldwprot.index).first()

repwprot.head()

oldwprot.head()

reprna.head()
len(nwprot)
len(nwprot[nwprot['GammaP']<.05])
sns.distplot(nwprot['GammaP'])

sns.distplot(wprot_res['GammaP'])
plt.savefig('wprot_res.pdf')

sns.distplot(cprot_res['GammaP'])
plt.savefig('cprot_res.pdf')


sns.distplot(bwprot[bwprot['GammaP']<.05]['Phase'])

sns.distplot(brna[brna['GammaP']<.05]['Phase'])

sns.distplot(brna['GammaP'])

len(brna)
len(brna[brna['GammaP']<.05])

sns.distplot(bcprot[bcprot['GammaP']<.05]['Phase'])

sns.distplot(bcprot['GammaP'])

len(bcprot)
len(bcprot[bcprot['GammaP']<.05])

a,b, _, _= scipy.stats.beta.fit(bwprot['GammaP'].values, floc=0, fscale=1)

plt.plot(np.linspace(0,1,1000),beta.pdf(np.linspace(0,1,1000),a,b))

a,b, _, _= scipy.stats.beta.fit(bcprot['GammaP'].values, floc=0, fscale=1)

plt.plot(np.linspace(0,1,1000),beta.pdf(np.linspace(0,1,1000),a,b))

sns.distplot(bwprot['GammaP'])

len(bwprot[bwprot['GammaP']>.5])
len(bwprot)
4747 - (1184*2)

len(brna[brna['GammaP']>.5])
len(brna)
9728 - (1979*2)
brna.ix['NCU03686']
brna.ix['NCU05154']

4736
sns.distplot(bwprot['GammaP'])

len(bwprot)
len(bwprot[bwprot['GammaP']<.05])


sns.jointplot(x='GammaP',y='Max_Amp',data=reprna)


best = repwprot[(repwprot['GammaP'] < .05) & (repwprot['Max_Amp'] > .5)]

len(best)

sns.jointplot(x='GammaP',y='Max_Amp',data=best)

sns.jointplot(x='GammaP',y='Max_Amp',data=repwprot[repwprot['GammaP']<.05])

sns.distplot(reprna['Phase'])

sns.distplot(reprna['GammaP'])

sns.distplot(repwprot['Phase'])

sns.distplot(repwprot['GammaP'])

sns.distplot(nwprot['GammaP'])

len(repwprot[repwprot['GammaP']<.05])

len(reprna[reprna['GammaP']<.05])
len(reprna[reprna['GammaP']>.95])

merged = pd.concat([bcprot['GammaP'], bwprot['GammaP']], axis=1,join='inner')
merged.columns = ['csp','wt']

sig = merged[(merged['csp'] < .05) & (merged['wt'] < .05)]

len(merged[merged['csp']<.05])
len(merged[merged['wt']<.05])

len(sig)

sigboth = merged[(merged['rna'] < .05) & (merged['prot'] < .05)]

len(sigboth)

sns.jointplot(x='csp',y='wt',data=merged)



sns.jointplot(x='csp',y='wt',data=merged,kind='kde',stat_func=None)
plt.savefig('cprotvswprot_pvals.pdf')

sns.jointplot(x='rna',y='prot',data=merged)


merged = pd.concat([brna['GammaP'], bwprot['GammaP']], axis=1,join='inner')
merged.columns = ['rna','prot']

sig = merged[(merged['rna'] < .05) & (merged['prot'] < .05)]

len(sig)

len(merged[merged['rna']<.05])

sigprot = merged[(merged['rna'] > .05) & (merged['prot'] < .05)]

len(sigprot)

sns.jointplot(x='rna',y='prot',data=merged,stat_func=None)
plt.savefig('rnavswprot_pvals_scatter.pdf')

len(merged[(merged['rna']>.5)&(merged['prot']<.05)])
len(merged)


sns.jointplot(x='rna',y='prot',data=merged,kind='kde',stat_func=None)
plt.savefig('rnavswprot_pvals.pdf')

sns.jointplot(x='rna',y='prot',data=merged)

merged = pd.concat([bcprot[bcprot['GammaP']<.05]['Phase'], bwprot[bwprot['GammaP']<.05]['Phase']], axis=1,join='inner')
merged.columns = ['csp','wt']

sns.jointplot(x='csp',y='wt',data=merged)

sns.jointplot(x='csp',y='wt',data=merged,kind='kde')

merged = pd.concat([nwprot['GammaP'], repwprot['GammaP']], axis=1,join='inner')
merged.columns = ['harm','rep']
merged.head()

sig = merged[(merged['harm'] < .05) | (merged['rep'] < .05)]

len(sig)

sigboth = merged[(merged['harm'] < .05) & (merged['rep'] < .05)]

len(sigboth)

sns.jointplot(x='harm',y='rep',data=sig)

len(nwprot[nwprot['GammaP']<.05])
len(nwprot)

merged = pd.concat([oldwprot['BH.Q'], repwprot['GammaP']], axis=1,join='inner')
merged.columns = ['old','new']

sig = merged[(merged['old'] < .05) | (merged['new'] < .05)]

len(sig)

sigboth = merged[(merged['old'] < .05) & (merged['new'] < .05)]

len(sigboth)

len(oldwprot[oldwprot['BH.Q']<.05])

sns.jointplot(x='old',y='new',data=sig)


merged = pd.concat([midwprot['GammaP'], repwprot['GammaP']], axis=1,join='inner')
merged.columns = ['old','new']

sig = merged[(merged['old'] < .05) | (merged['new'] < .01)]

len(sig)

sigboth = merged[(merged['old'] < .05) & (merged['new'] < .001)]

len(sigboth)

len(repwprot[repwprot['GammaP']<.001])

len(midwprot[midwprot['GammaP']<.05])

sns.jointplot(x='old',y='new',data=sig)


merged = pd.concat([circrepwprot['GammaP'], repwprot['GammaP']], axis=1,join='inner')
merged.columns = ['old','new']

sig = merged[(merged['old'] < .05) | (merged['new'] < .01)]

len(sig)

sigboth = merged[(merged['old'] < .01) & (merged['new'] < .001)]

len(sigboth)

len(repwprot[repwprot['GammaP']<.001])

len(circrepwprot[circrepwprot['GammaP']<.01])

sns.jointplot(x='old',y='new',data=sig)

merged = pd.concat([circrepsepwprot['GammaP'], repwprot['GammaP']], axis=1,join='inner')
merged.columns = ['old','new']

sns.jointplot(x='old',y='new',data=merged)

sig = merged[(merged['old'] < .05) | (merged['new'] < .05)]

len(sig)

sigboth = merged[(merged['old'] < .01) & (merged['new'] < .001)]

len(sigboth)

len(repwprot[repwprot['GammaP']<.001])

len(circrepsepwprot[circrepsepwprot['GammaP']<.01])

sns.jointplot(x='old',y='new',data=sig)

merged = pd.concat([circrepsepwprot['GammaP'], circrepwprot['GammaP']], axis=1,join='inner')
merged.columns = ['old','new']

sig = merged[(merged['old'] < .05) | (merged['new'] < .05)]

len(sig)

sigboth = merged[(merged['old'] < .05) & (merged['new'] < .05)]

len(sigboth)

len(repwprot[repwprot['GammaP']<.001])

len(circrepsepwprot[circrepsepwprot['GammaP']<.01])

sns.jointplot(x='old',y='new',data=sig)
sns.jointplot(x='old',y='new',data=merged)

bestnew = repwprot[(repwprot['GammaP'] < .5) & (repwprot['Max_Amp'] > .2)]

bestold = circrepsepwprot[(circrepsepwprot['GammaP'] < .5) & (circrepsepwprot['Max_Amp'] > .2)]


len(bestold)
len(bestnew)

merged = pd.concat([bestold['GammaP'], bestnew['GammaP']], axis=1,join='inner')
merged.columns = ['old','new']
sns.jointplot(x='old',y='new',data=merged


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

rna_old = pd.read_csv('JTK.wt_rna_scaled.txt',sep = '\t')
rna_old = rna_old.set_index('rownames(data)')


psandqs = pd.concat([ brna['Max_Amp'], bwprot['Max_Amp'],bcprot['Max_Amp'], rna_old['BH.Q'], rna_old['LAG'], brna['GammaP'], bwprot['GammaP'],bcprot['GammaP'],brna['Phase'], bwprot['Phase'], bcprot['Phase']], axis=1, join='inner')
psandqs.columns = ['rna_amp','wprot_amp','cprot_amp','rna_oldq', 'rna_old_phase', 'rna_pval','wprot_pval','cprot_pval','rna_phase','wprot_phase','cprot_phase']

sigboth = psandqs[(psandqs['wprot_pval']<.05) & (psandqs['cprot_pval'] < .05)]
len(psandqs[(psandqs['rna_pval']>.5) & (psandqs['wprot_pval'] < .05)])
len(psandqs)

brunner_genes = pd.read_csv('Brunner_genes.txt',sep='\t').columns.values.tolist()[0].split(' ')

set(brunner_genes).intersection(difs[(difs > 6) & (difs < 17)].index.tolist())
set(brunner_genes).intersection(difs[(difs < -6) & (difs > -17)].index.tolist())

len(set(brunner_genes).intersection(difs[(difs < -6) & (difs > -17)].index.tolist()))
len(set(brunner_genes).intersection(difs[(difs > 6) & (difs < 17)].index.tolist()))
len(set(brunner_genes).intersection(sigboth.index.tolist()))
len(sigboth)
len(brunner_genes)

#1.331257926 times as many as expected at random a lot of the csp-1 targets are not circadian though

.389830508

.292828685

difs = sigboth['wprot_phase'] - sigboth['cprot_phase']
len(difs[(difs > 6) & (difs < 17)].index)

len(difs[(difs < -6) & (difs > -17)].index)

difs[(difs > 6) & (difs < 17)].index.tolist() + difs[(difs < -6) & (difs > -17)].index.tolist()

theta = np.linspace(0,2*np.pi,145)
theta_matrix = np.meshgrid(theta)[0].tolist()
transformed = [(i % 23)*np.pi*2/23 for i in difs.values.tolist()]

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
plt.savefig('wtvscsp.pdf')

sigboth = psandqs[(psandqs['wprot_pval']<.05) & (psandqs['rna_pval'] < .05)]
len(sigboth)
len(psandqs)

difs = sigboth['wprot_phase'] - sigboth['rna_phase']
difs.values

theta = np.linspace(0,2*np.pi,145)
theta_matrix = np.meshgrid(theta)[0].tolist()
transformed = [(i % 23)*np.pi*2/23 for i in difs.values.tolist()]

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
plt.savefig('rnavsprot.pdf')


theta = np.linspace(0,2*np.pi,145)
theta_matrix = np.meshgrid(theta)[0].tolist()
transformed = [(i % 23)*np.pi*2/23 for i in psandqs[psandqs['wprot_pval']<.05]['wprot_phase'].values.tolist()]

transformed = transformed + [(i - (2*np.pi)) for i in transformed] + [(i + (2*np.pi)) for i in transformed]

kernel = stats.gaussian_kde(transformed + [(i + (2*np.pi)) for i in transformed])
Z = kernel(theta_matrix)

kde = KernelDensity(kernel='gaussian',bandwidth=0.25).fit(np.asarray(transformed).reshape(-1,1))
Z = kde.score_samples(np.asarray(theta_matrix).reshape(-1,1))

Z  = scale(Z)
Z = Z - np.min(Z)
min_max_scaler = MinMaxScaler()
Z = min_max_scaler.fit_transform(Z.reshape(-1,1))

transformedn = [(i % 23)*np.pi*2/23 for i in psandqs[psandqs['cprot_pval']<.05]['cprot_phase'].values.tolist()]
transformedn = transformedn + [(i - (2*np.pi)) for i in transformedn] + [(i + (2*np.pi)) for i in transformedn]

kden = KernelDensity(kernel='gaussian',bandwidth=0.25).fit(np.asarray(transformedn).reshape(-1,1))
Zn = kden.score_samples(np.asarray(theta_matrix).reshape(-1,1))


Zn = scale(Zn)
Zn = Zn - np.min(Zn)
min_max_scaler = MinMaxScaler()
Zn = min_max_scaler.fit_transform(Zn.reshape(-1,1))

transformedm = [(i % 23)*np.pi*2/23 for i in psandqs[psandqs['rna_pval']<.05]['rna_phase'].values.tolist()]
transformedm = transformedm + [(i - (2*np.pi)) for i in transformedm] + [(i + (2*np.pi)) for i in transformedm]

kdem = KernelDensity(kernel='gaussian',bandwidth=0.25).fit(np.asarray(transformedm).reshape(-1,1))
Zm = kdem.score_samples(np.asarray(theta_matrix).reshape(-1,1))


Zm = scale(Zm)
Zm = Zm - np.min(Zm)
min_max_scaler = MinMaxScaler()
Zm = min_max_scaler.fit_transform(Zm.reshape(-1,1))

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(260)
CS = plt.plot(np.linspace(0,2*np.pi,145), Z, c= 'b', label='WT Protein')
CS = plt.plot(np.linspace(0,2*np.pi,145), Zn, c= 'r', label='Delta CSP-1 Protein')
CS = plt.plot(np.linspace(0,2*np.pi,145), Zm, c= 'g', label='WT RNA')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('densitieseach.pdf')

exp(0.3143)

sigall = psandqs[(psandqs['wprot_pval']<.05) & (psandqs['cprot_pval'] < .05) & (psandqs['rna_pval'] < .05)]

theta = np.linspace(0,2*np.pi,145)
theta_matrix = np.meshgrid(theta)[0].tolist()
transformed = [(i % 23)*np.pi*2/23 for i in sigall['wprot_phase'].values.tolist()]

transformed = transformed + [(i - (2*np.pi)) for i in transformed] + [(i + (2*np.pi)) for i in transformed]

kernel = stats.gaussian_kde(transformed + [(i + (2*np.pi)) for i in transformed])
Z = kernel(theta_matrix)

kde = KernelDensity(kernel='gaussian',bandwidth=0.25).fit(np.asarray(transformed).reshape(-1,1))
Z = kde.score_samples(np.asarray(theta_matrix).reshape(-1,1))

Z  = scale(Z)
Z = Z - np.min(Z)
min_max_scaler = MinMaxScaler()
Z = min_max_scaler.fit_transform(Z.reshape(-1,1))

transformedn = [(i % 23)*np.pi*2/23 for i in sigall['cprot_phase'].values.tolist()]
transformedn = transformedn + [(i - (2*np.pi)) for i in transformedn] + [(i + (2*np.pi)) for i in transformedn]

kden = KernelDensity(kernel='gaussian',bandwidth=0.25).fit(np.asarray(transformedn).reshape(-1,1))
Zn = kden.score_samples(np.asarray(theta_matrix).reshape(-1,1))


Zn = scale(Zn)
Zn = Zn - np.min(Zn)
min_max_scaler = MinMaxScaler()
Zn = min_max_scaler.fit_transform(Zn.reshape(-1,1))

transformedm = [(i % 23)*np.pi*2/23 for i in sigall['rna_phase'].values.tolist()]
transformedm = transformedm + [(i - (2*np.pi)) for i in transformedm] + [(i + (2*np.pi)) for i in transformedm]

kdem = KernelDensity(kernel='gaussian',bandwidth=0.25).fit(np.asarray(transformedm).reshape(-1,1))
Zm = kdem.score_samples(np.asarray(theta_matrix).reshape(-1,1))


Zm = scale(Zm)
Zm = Zm - np.min(Zm)
min_max_scaler = MinMaxScaler()
Zm = min_max_scaler.fit_transform(Zm.reshape(-1,1))

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(115)
CS = plt.plot(np.linspace(0,2*np.pi,145), Z, c= 'b', label='WT Protein')
CS = plt.plot(np.linspace(0,2*np.pi,145), Zn, c= 'r', label='Delta CSP-1 Protein')
CS = plt.plot(np.linspace(0,2*np.pi,145), Zm, c= 'g', label='WT RNA')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('densitiesall.pdf')




forKDE = pol2cart(psandqs[psandqs['rna_pval']<.05]['rna_oldq'].values,psandqs[psandqs['rna_pval']<.05]['rna_phase'].values*np.pi*2/23)

values = np.vstack([forKDE[0], forKDE[1]])

R = 1
r = np.linspace(0,R,100)
theta = np.linspace(-np.pi,np.pi,200)
radius_matrix, theta_matrix = np.meshgrid(r,theta)
positions = np.vstack([radius_matrix.ravel(), theta_matrix.ravel()])

cart_pos = pol2cart(positions[0],positions[1])

cart_pos = np.array([cart_pos[0],cart_pos[1]])

kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(cart_pos).T, radius_matrix.shape)

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(70)
CS = plt.contour(theta_matrix, radius_matrix, Z,
                cmap=cm.plasma)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('Old RNA Q-value vs New RNA Phase')
plt.savefig('RNA_new_phase.pdf')


forKDE = pol2cart(psandqs[psandqs['wprot_pval']<.05]['rna_pval'].values,psandqs[psandqs['wprot_pval']<.05]['wprot_phase'].values*np.pi*2/23)

values = np.vstack([forKDE[0], forKDE[1]])

R = 1
r = np.linspace(0,R,100)
theta = np.linspace(-np.pi,np.pi,200)
radius_matrix, theta_matrix = np.meshgrid(r,theta)
positions = np.vstack([radius_matrix.ravel(), theta_matrix.ravel()])

cart_pos = pol2cart(positions[0],positions[1])

cart_pos = np.array([cart_pos[0],cart_pos[1]])

kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(cart_pos).T, radius_matrix.shape)

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(70)
CS = plt.contour(theta_matrix, radius_matrix, Z,
                cmap=cm.plasma)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('WTprot Pval vs RNA Phase')
plt.savefig('wprotpval_v_rnaphase.pdf')

forKDE = pol2cart(psandqs[psandqs['cprot_pval']<.05]['wprot_pval'].values,psandqs[psandqs['cprot_pval']<.05]['cprot_phase'].values*np.pi*2/23)

values = np.vstack([forKDE[0], forKDE[1]])

R = 1
r = np.linspace(0,R,100)
theta = np.linspace(-np.pi,np.pi,200)
radius_matrix, theta_matrix = np.meshgrid(r,theta)
positions = np.vstack([radius_matrix.ravel(), theta_matrix.ravel()])

cart_pos = pol2cart(positions[0],positions[1])

cart_pos = np.array([cart_pos[0],cart_pos[1]])

kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(cart_pos).T, radius_matrix.shape)

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(70)
CS = plt.contour(theta_matrix, radius_matrix, Z,
                cmap=cm.plasma)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('WTprot Pval vs Delta CSP-1 Phase')
plt.savefig('wprotpval_v_cspphase.pdf')

forKDE = pol2cart(psandqs[psandqs['wprot_pval']<.05]['wprot_pval'].values,psandqs[psandqs['wprot_pval']<.05]['wprot_phase'].values*np.pi*2/23)

values = np.vstack([forKDE[0], forKDE[1]])

R = 0.05
r = np.linspace(0,R,100)
theta = np.linspace(-np.pi,np.pi,200)
radius_matrix, theta_matrix = np.meshgrid(r,theta)
positions = np.vstack([radius_matrix.ravel(), theta_matrix.ravel()])

cart_pos = pol2cart(positions[0],positions[1])

cart_pos = np.array([cart_pos[0],cart_pos[1]])

kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(cart_pos).T, radius_matrix.shape)

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(70)
CS = plt.contour(theta_matrix, radius_matrix, Z,
                cmap=cm.plasma)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('WTprot Pval vs WTprot Phase')
plt.savefig('wprotpval_v_wprotphase.pdf')

min_max_scaler = MinMaxScaler()
scaled = [i[0] for i in min_max_scaler.fit_transform(psandqs[psandqs['cprot_pval']<.05]['cprot_amp'].reshape(-1,1))]

forKDE = pol2cart(scaled,psandqs[psandqs['cprot_pval']<.05]['cprot_phase'].values*np.pi*2/23)

values = np.vstack([forKDE[0], forKDE[1]])

np.shape(values)

R = 1
r = np.linspace(0,R,100)
theta = np.linspace(-np.pi,np.pi,200)
radius_matrix, theta_matrix = np.meshgrid(r,theta)
positions = np.vstack([radius_matrix.ravel(), theta_matrix.ravel()])

cart_pos = pol2cart(positions[0],positions[1])

cart_pos = np.array([cart_pos[0],cart_pos[1]])

kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(cart_pos).T, radius_matrix.shape)

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(70)
CS = plt.contour(theta_matrix, radius_matrix, Z,
                cmap=cm.plasma)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('CSPprot Amp vs CSPprot Phase')
plt.savefig('cspamp_v_cprotphase.pdf')


min_max_scaler = MinMaxScaler()
scaled = [i[0] for i in min_max_scaler.fit_transform(psandqs[psandqs['wprot_pval']<.05]['wprot_amp'].reshape(-1,1))]

forKDE = pol2cart(scaled,psandqs[psandqs['wprot_pval']<.05]['wprot_phase'].values*np.pi*2/23)

values = np.vstack([forKDE[0], forKDE[1]])

R = 1
r = np.linspace(0,R,100)
theta = np.linspace(-np.pi,np.pi,200)
radius_matrix, theta_matrix = np.meshgrid(r,theta)
positions = np.vstack([radius_matrix.ravel(), theta_matrix.ravel()])

cart_pos = pol2cart(positions[0],positions[1])

cart_pos = np.array([cart_pos[0],cart_pos[1]])

kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(cart_pos).T, radius_matrix.shape)

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(70)
CS = plt.contour(theta_matrix, radius_matrix, Z,
                cmap=cm.plasma)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('WTprot Amp vs WTprot Phase')
plt.savefig('wtamp_v_wprotphase.pdf')


min_max_scaler = MinMaxScaler()
scaled = [i[0] for i in min_max_scaler.fit_transform(psandqs[psandqs['rna_pval']<.05]['rna_amp'].reshape(-1,1))]

forKDE = pol2cart(scaled,psandqs[psandqs['rna_pval']<.05]['rna_phase'].values*np.pi*2/23)

values = np.vstack([forKDE[0], forKDE[1]])

R = 1
r = np.linspace(0,R,100)
theta = np.linspace(-np.pi,np.pi,200)
radius_matrix, theta_matrix = np.meshgrid(r,theta)
positions = np.vstack([radius_matrix.ravel(), theta_matrix.ravel()])

cart_pos = pol2cart(positions[0],positions[1])

cart_pos = np.array([cart_pos[0],cart_pos[1]])

kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(cart_pos).T, radius_matrix.shape)

plt.figure()
ax = plt.axes(polar=True)
xtickslbls = ['0', '2.875','5.75','8.625','11.5','14.375','17.25','20.125']
ax.set_xticklabels(xtickslbls)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2.0)
ax.set_rlabel_position(70)
CS = plt.contour(theta_matrix, radius_matrix, Z,
                cmap=cm.plasma)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('RNA Amp vs RNA Phase')
plt.savefig('rnaamp_v_rnaphase.pdf')


prot_sig = psandqs[psandqs['wprot_pval']<.05]
prot_sig = prot_sig[prot_sig['cprot_pval']<.05]
prot_sig['pdelt'] = prot_sig['wprot_phase'] - prot_sig['cprot_phase']
prot_sig['pdelt1'] = np.sin(prot_sig['pdelt']*np.pi*2/23)
prot_sig['pdelt2'] = np.cos(prot_sig['pdelt']*np.pi*2/23)

graph = sns.jointplot(prot_sig.pdelt1, prot_sig.pdelt2,stat_func=None,color="r", shade_lowest=False,space=0,kind = 'kde')
plt.savefig('jen_wt_v_csp_phase.pdf')

sig =[]
sig = psandqs[psandqs['wprot_pval']<.05]
sig = sig[sig['rna_pval']<.05]
sig['pdelt'] = sig['wprot_phase'] - sig['rna_phase']
sig['pdelt1'] = np.sin(sig['pdelt']*np.pi*2/23)
sig['pdelt2'] = np.cos(sig['pdelt']*np.pi*2/23)

graph = sns.jointplot(sig.pdelt1, sig.pdelt2,kind = 'kde',color="g", shade_lowest=False,space=0,stat_func=None)
plt.savefig('jen_prot_v_rna_phase.pdf')


w_sig=[]
w_sig = psandqs[psandqs['wprot_pval']<.05]
w_sig['wt1'] = np.sin(w_sig['wprot_phase']*np.pi*2/23)
w_sig['wt2'] = np.cos(w_sig['wprot_phase']*np.pi*2/23)


graph = sns.jointplot(w_sig.wt1, w_sig.wt2,kind = 'kde',color="b", shade_lowest=False,space=0,stat_func=None)
plt.savefig('jen_wprot.pdf')

c_sig=[]
c_sig = psandqs[psandqs['cprot_pval']<.05]
c_sig['ct1'] = np.sin(c_sig['cprot_phase']*np.pi*2/23)
c_sig['ct2'] = np.cos(c_sig['cprot_phase']*np.pi*2/23)

graph = sns.jointplot(c_sig.ct1, c_sig.ct2,kind = 'kde',color="b", shade_lowest=False,space=0,stat_func=None)
plt.savefig('jen_cprot.pdf')

r_sig=[]
r_sig = psandqs[psandqs['rna_pval']<.05]
r_sig['rt1'] = np.sin(r_sig['rna_phase']*np.pi*2/23)
r_sig['rt2'] = np.cos(r_sig['rna_phase']*np.pi*2/23)

graph = sns.jointplot(r_sig.rt1, r_sig.rt2,kind = 'kde',color="b", shade_lowest=False,space=0,stat_func=None)
plt.savefig('jen_rna.pdf')
