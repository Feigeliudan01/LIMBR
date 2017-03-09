import numpy as np, pandas as pd, matplotlib.pyplot as plt; plt.rc('text', usetex=True)
import seaborn as sns; sns.set(style="white", context="talk")
import scipy.stats as st
from sklearn.metrics import mutual_info_score
from  matplotlib import colors, ticker, cm
import scipy.stats as stats
from sklearn import mixture
from astroML.utils import convert_2D_cov
from astroML.plotting.tools import draw_ellipse
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=True)
from matplotlib.colors import BoundaryNorm
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
circularlib=importr('circular')
circular=robjects.r('circular')
corcircular=robjects.r('cor.circular')
from sklearn.cluster import KMeans
from sklearn import preprocessing
from sklearn.decomposition import PCA
import scipy.cluster.hierarchy as sch
from palettable.colorbrewer.diverging import RdBu_4

wt = pd.read_csv('wt_lowess_normed__jtkout_GammaP.txt',sep='\t',index_col=0)
csp = pd.read_csv('csp_lowess_normed__jtkout_GammaP.txt',sep='\t',index_col=0)
rna = pd.read_csv('rna_lowess_normed__jtkout_GammaP.txt',sep='\t',index_col=0)

prot = pd.merge(wt[['GammaP','Phase']], csp[['GammaP','Phase']], right_index=True, left_index=True,suffixes=('_wt','_csp'))

prot['zscore_wt'] = st.norm.ppf(1-prot['GammaP_wt'])
prot['zscore_csp'] = st.norm.ppf(1-prot['GammaP_csp'])
prot['delta_z'] = prot['zscore_wt'] - prot['zscore_csp']
prot['delta_p'] = prot['Phase_csp'] - prot['Phase_wt']


rp = pd.merge(wt[['GammaP','Phase']], rna[['GammaP','Phase']], right_index=True, left_index=True,suffixes=('_wt','_rna'))

rp['zscore_wt'] = st.norm.ppf(1-rp['GammaP_wt'])
rp['zscore_rna'] = st.norm.ppf(1-rp['GammaP_rna'])
rp['delta_p'] = rp['Phase_wt'] - rp['Phase_rna']

tau, p_value = stats.kendalltau(prot['zscore_wt'].values, prot['zscore_csp'].values)
r, p_value = stats.pearsonr(prot['zscore_wt'].values, prot['zscore_csp'].values)


###########Hexbin Plot of Z score Distribution

ax = plt.subplot(111)
hb = ax.hexbin(prot['zscore_wt'].values, prot['zscore_csp'].values, gridsize=48, extent=(-2,4,-2,4))
ax.set_xlabel(r"WT Z Score", fontsize=12)
ax.set_ylabel(r'$\Delta$CSP-1 Z Score', fontsize=12)
cb = fig.colorbar(hb, ax=ax, ticks=np.linspace(0,20,11))
cb.set_label('count', fontsize=12)
ax.set_title('Distribution of Z Scores from eJTK of Proteomics Data by Genotype', fontsize=16)
ax.text(0.015, 0.985, "Pearson's r = %f, p = %e"%(r, p_value), transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=dict(facecolor='white'))
plt.savefig('prot_hex.pdf')
plt.close()

norm_map = gen_norm_dict(wt.columns.values)
wt_pnorm = pool_normalize(wt,norm_map)
wt_pnorm = wt_pnorm.replace([np.inf, -np.inf], np.nan)
wt_pnorm = wt_pnorm.dropna()
wt_pnorm = wt_pnorm.sort_index(axis=1)
wt_pnorm = qnorm(wt_pnorm)
wt_pnorm = wt_pnorm.groupby(level='Protein').mean()
wt_pnorm.index.names = ['#']
wt_pnorm.to_csv('wt_qnorm.txt',sep='\t')

norm_map = gen_norm_dict(csp.columns.values)
csp_pnorm = pool_normalize(csp,norm_map)
csp_pnorm = csp_pnorm.replace([np.inf, -np.inf], np.nan)
csp_pnorm = csp_pnorm.dropna()
csp_pnorm = csp_pnorm.sort_index(axis=1)
csp_pnorm = qnorm(csp_pnorm)
csp_pnorm = csp_pnorm.groupby(level='Protein').mean()
csp_pnorm.index.names = ['#']
csp_pnorm.to_csv('csp_qnorm.txt',sep='\t')

wt_unnorm = pd.read_csv('wt_qnorm__jtkout_GammaP.txt',sep='\t',index_col=0)
csp_unnorm = pd.read_csv('csp_qnorm__jtkout_GammaP.txt',sep='\t',index_col=0)
unnorm = pd.merge(wt_unnorm[['GammaP','Phase']], csp_unnorm[['GammaP','Phase']], right_index=True, left_index=True,suffixes=('_wt','_csp'))

unnorm['zscore_wt'] = st.norm.ppf(1-unnorm['GammaP_wt'])
unnorm['zscore_csp'] = st.norm.ppf(1-unnorm['GammaP_csp'])
unnorm['delta_z'] = unnorm['zscore_wt'] - unnorm['zscore_csp']
unnorm['delta_p'] = unnorm['Phase_csp'] - unnorm['Phase_wt']

r, p_value = stats.pearsonr(unnorm['zscore_wt'].values, unnorm['zscore_csp'].values)
tau, p_value = stats.kendalltau(unnorm['zscore_wt'].values, unnorm['zscore_csp'].values)

ax = plt.subplot(111)
hb = ax.hexbin(unnorm['zscore_wt'].values, unnorm['zscore_csp'].values, gridsize=48, extent=(-2,4,-2,4))
ax.set_xlabel(r"WT Z Score", fontsize=12)
ax.set_ylabel(r'$\Delta$CSP-1 Z Score', fontsize=12)
cb = fig.colorbar(hb, ax=ax, ticks=np.linspace(0,20,11))
cb.set_label('count', fontsize=12)
ax.set_title('Distribution of Z Scores from eJTK of Proteomics Data by Genotype', fontsize=16)
ax.text(0.015, 0.985, "Pearson's r = %f, p = %e"%(r, p_value), transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=dict(facecolor='white'))
plt.savefig('unnorm_hex.pdf')
plt.close()


X = np.asarray([prot['zscore_wt'].values, prot['zscore_csp'].values]).T
Nclusters = np.arange(1, 8)
clfs = [mixture.GaussianMixture(n_components=N, covariance_type='full').fit(X) for N in Nclusters]
BICs = np.array([clf.bic(X) for clf in clfs])
print("convergence:", [clf.converged_ for clf in clfs])

# plot the BIC
ax = plt.subplot(111)
ax.plot(Nclusters, BICs / len(X), '-', c='k')
ax.set_xlabel(r'$n. clusters$', fontsize=12)
ax.set_ylabel(r'$BIC / n. observations$', fontsize=12)
plt.title(r'Determination of \# of Components From BIC', fontsize=16)
plt.savefig('BIC.pdf')
plt.close()


clf = clfs[np.argmin(BICs)]
#log_dens = clf.score(Xgrid).reshape((70, 70))

###########Gaussian Mixture Model

# scatter the points
ax = plt.subplot(111)
ax.plot(X[:, 0], X[:, 1], ',k', alpha=0.3, zorder=1)
for i in range(clf.n_components):
    mean = clf.means_[i]
    cov = clf.covariances_[i]
    if cov.ndim == 1:
        cov = np.diag(cov)
    draw_ellipse(mean, cov, ax=ax, fc='none', ec=['b','r'][i], zorder=2)

plt.plot([], label='Predicted Distribution of Circadian Proteins', color="blue",linewidth=.5)
plt.plot([], label='Predicted Distribution of Non-Circadian Proteins', color="red",linewidth=.5)
plt.plot([], label='Probability of Being Drawn From Circadian Distribution', color="black",linewidth=.5)
plt.plot([], label='Z Score Corresponding to p=0.05', color="black",linewidth=.5, linestyle='--')
plt.legend(loc="lower left")
ax.set_xlabel(r"WT Z Score", fontsize=12)
ax.set_ylabel(r'$\Delta$CSP-1 Z Score', fontsize=12)
plt.title('Gaussian Mixture Model of Z Score with 2 Components', fontsize=16)
ax.scatter(X[:, 0], X[:, 1],color='k', alpha=.3, zorder=1,s=1.5)
x, y = np.meshgrid(np.linspace(-4., 5.,1000),np.linspace(-4., 5.,1000))
XX = np.array([x.ravel(), y.ravel()]).T
Z = clf.predict_proba(XX)[:,0]
Z = Z.reshape(x.shape)
CS = ax.contour(x, y, Z, color='k', levels=np.linspace(.9, 1, 11),linewidths=.75)
#manual_locations = list(zip(np.linspace(-2,-1,8),np.linspace(4,4.7,8)))
ax.axvline(x=st.norm.ppf(.95), linestyle='--', color='black', linewidth=0.4, zorder=3)
ax.axhline(y=st.norm.ppf(.95), linestyle='--', color='black', linewidth=0.4, zorder=3)
manual_locations = [(3.8,-2),(4,-2),(4.1,-2),(4.2,-2),(4.4,-2),(4.6,-2),(4.9,-2)]
plt.clabel(CS, manual=manual_locations)
plt.xlim(-4,5)
plt.ylim(-4,5)
plt.savefig('GMM.pdf')
plt.close()




defcirc = prot[clf.predict_proba(X)[:,0] > .95]
#defcirc = prot[(prot['GammaP_wt']<.05)&(prot['GammaP_csp']<.05)]

###########Phase Lag

cpw = circular(robjects.FloatVector([j*2*np.pi/22 for j in defcirc['Phase_wt'].values]),units="radians", zero=0, rotation='clock')
cpc = circular(robjects.FloatVector([j*2*np.pi/22 for j in defcirc['Phase_csp'].values]),units="radians", zero=0, rotation='clock')

circorr = corcircular(cpw,cpc,test='TRUE')

circorr[0].r_repr()

def phase_hist(phases,ptitle):
    N = 11
    phases = [i%22 for i in phases]
    radii = np.histogram(phases,N)[0]
    bottom = 0
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    width = (2*np.pi) / N
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(0.5*np.pi)
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(180)
    ax.set_xticks(np.pi/180. * np.linspace(0,  360, N, endpoint=False))
    ax.set_xticklabels([str(i)+' hrs' for i in range(0,2*N+1,2)])
    bars = ax.bar(theta, radii, width=width, bottom=bottom)
    # Use custom colors and opacity
    for r, bar in zip(radii, bars):
        bar.set_facecolor('black')
        bar.set_alpha(0.8)
    plt.title(ptitle, y=1.08, fontsize=16)
    ax.text(0.85, 0.05, r'\rho_{cc} = %f'"\n"r'p = %e'%(float(circorr[0].r_repr()), float(circorr[2].r_repr())), transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=dict(facecolor='white'))
    plt.tight_layout()
    plt.savefig('WT_CSP_phasehist.pdf')
    plt.close()


phase_hist(defcirc['delta_p'].values,r'Histogram of Phase Lags Between'"\n"r'Expression in WT and $\Delta$CSP-1 Genotypes'"\n"r'For Proteins Circadian in Both Genotypes')

#################


#need to write function which does tsplot for WT and delta csp
#then need to make these plots for the highest probability 10 Proteins
#then need to make these plots for the 10 proteins with the highest delta z score

wt_exp = pd.read_csv('wt_lowess_normed.txt',sep='\t',index_col=0)
csp_exp = pd.read_csv('csp_lowess_normed.txt',sep='\t',index_col=0)
wt_exp_circ = wt_exp.loc[prot[prot['GammaP_wt']<.05].index.values]
csp_exp_circ = csp_exp.loc[prot[prot['GammaP_csp']<.05].index.values]

def scale_df(df):
    sc = df.apply(lambda x: x-np.mean(x),axis=1)
    sc = sc.apply(lambda x: x/np.std(x),axis=1)
    return sc

scaled_wt = scale_df(wt_exp_circ)
scaled_csp = scale_df(csp_exp_circ)

z_min, z_max = -np.abs(scaled_wt).max(), np.abs(scaled_wt).max()

random_state = 42
y_pred = KMeans(n_clusters=2, random_state=random_state).fit_predict(scaled_wt)


def plot_heatmap(df,fname,ptitle):
    z_min, z_max = -np.abs(df).max(), np.abs(df).max()
    random_state = 42
    y_pred = KMeans(n_clusters=2, random_state=random_state).fit_predict(df)
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ordered = reorder_heatmap(df)
    pc = ax.pcolor(ordered, cmap='seismic', vmin=z_min*.8, vmax=z_max*.8)
    for i in range(0,ordered.shape[1],3):
        ax.axvline(x=i, linestyle='-', color='black', linewidth=0.4,    zorder=3)
    ax.set_xticks([i+1.5 for i in range(0,ordered.shape[1],3)])
    ax.set_xticklabels([i*2 for i in range(1,(int(ordered.shape[1]/3)+1))])
    ax.set_yticks([])
    ax.set_title(ptitle, fontsize=16)
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
    cb = fig.colorbar(pc, cax=cbar_ax)
    cbar_ax.text(2.5,.49,r'$\sigma$', fontsize=12)
    plt.savefig(fname+'.pdf')

plot_heatmap(scaled_wt.values,r'WT_circ_genes','Expression as a Function of Timepoint in Circadian Genes in the WT Strain')
plot_heatmap(scaled_csp.values,r'CSP_circ_genes',' Expression as a Function of Timepoint in Circadian Genes in the $\Delta$CSP-1 Strain')

def reorder_heatmap(arr):
    Y = sch.linkage(arr, method='complete', metric='euclidean')
    Z = sch.dendrogram(Y, orientation='left',no_plot=True)
    idx = Z['leaves']
    return arr[idx,:]

plot_heatmap(scaled_wt.values)

def gen_tsplot(df,fname):
    ldf = pd.melt(df.reset_index(), id_vars=[df.reset_index().columns.values.tolist()[0]], value_vars=df.columns.values.tolist())
    ldf['replicate'] = ldf['variable'].apply(lambda x: x.split('.')[1] if len(x.split('.')) != 1 else '0')
    ldf['variable'] = ldf['variable'].apply(lambda x: x.split('.')[0][2:])
    ldf.columns = ['Gene','Timepoint','Expression','Replicate']
    sns.tsplot(time="Timepoint", value="Expression", unit='Replicate', condition="Gene",data=ldf,color="Set1",ci=95)
    plt.savefig(fname+'.pdf')
    plt.close()


gen_tsplot(scaled_wt.loc[['NCU08312','NCU04238']],'test')

def plot_cormat(df,fname,ptitle):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    cors = np.corrcoef(df.values.T)
    pc = ax.pcolor(cors, cmap='seismic')
    for i in range(0,cors.shape[1],3):
        ax.axvline(x=i, linestyle='-', color='black', linewidth=0.7,    zorder=3)
    for i in range(0,cors.shape[0],3):
        ax.axhline(y=i, linestyle='-', color='black', linewidth=0.7,    zorder=3)
    ax.set_xticks([i+1.5 for i in range(0,cors.shape[1],3)])
    ax.set_xticklabels([i*2 for i in range(1,(int(cors.shape[1]/3)+1))])
    ax.set_yticks([i+1.5 for i in range(0,cors.shape[1],3)])
    ax.set_yticklabels([i*2 for i in range(1,(int(cors.shape[1]/3)+1))])
    ax.set_title(ptitle, fontsize=16)
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
    fig.colorbar(pc, cax=cbar_ax)
    plt.savefig(fname+'.pdf')
    plt.close()


plot_cormat(scaled_wt,'final_cormat','Correlation Matrix of LIMBR Normalized WT Samples Ordered by Timepoint')

wt_qnorm = pd.read_csv('wt_qnorm.txt',sep='\t',index_col=0)

plot_cormat(wt_qnorm,'initial_cormat','Correlation Matrix of Traditionally Normalized WT Samples Ordered by Timepoint')

m_x

pca = PCA(n_components=4)
X_r = pca.fit(scaled_wt.values.T).transform(scaled_wt.values.T)

plt.scatter(X_r[:,0],X_r[:,1])
plt.show()





def plt_line(arr,ind_1,ind_2):
    plt.plot([arr[ind_1][2],arr[ind_2][2]],[arr[ind_1][3],arr[ind_2][3]],color='k')


for i in range(int(X_r.shape[0]/3)):
    plt_line(X_r,(i*3),(i*3+1))
    plt_line(X_r,(i*3),(i*3+2))
    plt_line(X_r,(i*3+1),(i*3+2))

plt.show()










wt_class = pd.read_csv('wt_classes.csv')
csp_class = pd.read_csv('csp_classes.csv')

wt_trends = pd.read_csv('wt_lowess_normed_trends.txt',sep='\t')
wt_trends = wt_trends.drop('Unnamed: 0',axis=1)
wt_trends['Trend'] = range(1,len(wt_trends)+1)
wt_trends = pd.melt(wt_trends, id_vars=['Trend'], value_vars=wt_trends.columns.values.tolist()[:-1])

csp_trends = pd.read_csv('csp_lowess_normed_trends.txt',sep='\t')
csp_trends = csp_trends.drop('Unnamed: 0',axis=1)
csp_trends['Trend'] = range(1,len(csp_trends)+1)
csp_trends = pd.melt(csp_trends, id_vars=['Trend'], value_vars=csp_trends.columns.values.tolist()[:-1])


def get_fstats(ts,cs,ctypes):
    results = {}
    for class_val in ctypes:
        groups = []
        for i in list(set(cs[class_val].values)):
            groups.append(list(cs[cs[class_val]==i]['Sample ID'].values))
        fortest = {}
        for j in np.sort(list(set(ts['Trend'].values))):
            templist = []
            for i in range(len(groups)):
                templist.append(list(ts[ts['variable'].isin(groups[i]) & (ts['Trend'] == j)]['value'].values))
            fortest[j] = templist
        tempdict = {}
        for k, v in fortest.items():
            tempdict[k] = stats.f_oneway(*v)[1]
        results[class_val] = tempdict
    return pd.DataFrame.from_dict(results)


wt_bias = get_fstats(wt_trends,wt_class,['TMT Set', 'HpH Fractionation','Post-Digest SPE', 'Post-Digest BCA', 'Digest'])
csp_bias = get_fstats(csp_trends,csp_class,['TMT Set', 'HpH Fractionation','Post-Digest SPE', 'Post-Digest BCA', 'Digest'])


def gen_biasmap(df):
    plt.pcolor(df.values)
    for y in range(df.values.shape[0]):
        for x in range(df.values.shape[1]):
            if df.values[y,x] < .05/df.values.shape[1]:
                plt.text(x + 0.5, y + 0.5, '*', horizontalalignment='center',verticalalignment='center')
    plt.show()

gen_biasmap(wt_bias)



#plt.colorbar()
plt.show()


plt.pcolor(scaled_wt[y_pred==0], cmap='seismic', vmin=z_min, vmax=z_max)
plt.colorbar()
plt.show()


scaled_csp = preprocessing.scale(csp_exp_circ.values, axis=1)
z_min, z_max = -np.abs(scaled_csp).max(), np.abs(scaled_csp).max()


plt.pcolor(scaled_csp[y_pred==1], cmap='seismic', vmin=z_min, vmax=z_max)
plt.colorbar()
plt.show()

plt.pcolor(scaled_csp[y_pred==0], cmap='seismic', vmin=z_min, vmax=z_max)
plt.colorbar()
plt.show()










for i in range(clf.n_components):
    Z = clf.predict_proba(XX)[:,i]
    Z = Z.reshape(x.shape)
    CS.append(ax.contour(x, y, Z, colors=['r','b'][i], levels=np.linspace(0, 1, 100)))

plt.show()

plt.plot([], label='Circadian Proteins', color="blue")
plt.plot([], label='Non-Circadian Proteins', color="red")
plt.legend(loc="upper left")
ax.set_xlabel(r"WT Z Score")
ax.set_ylabel(r'$\Delta$ CSP Z Score')
plt.title('Gaussian Mixture Model of Z Score with 2 Components')
plt.show()
# plot the components


# label the plot
ax.text(0.05, 0.95, "N = %i points" % Npts,ha='left', va='top', transform=ax.transAxes, bbox=dict(fc='w', ec='k'))

    ax.set_xlim(-5, 105)
    ax.set_ylim(-5, 105)


ax_list[0].xaxis.set_major_formatter(plt.NullFormatter())
ax_list[2].yaxis.set_major_formatter(plt.NullFormatter())

for i in (0, 1):
    ax_list[i].set_ylabel('$y$')

for j in (1, 2):
    ax_list[j].set_xlabel('$x$')

ax_list[-1].legend(loc=1)

ax_list[-1].set_xlabel('n. clusters')
ax_list[-1].set_ylabel('$BIC / N$')
ax_list[-1].set_ylim(16, 18.5)

plt.show()









def phase_hist(df,fname):
    phases = df['delta_p'].values.tolist()
    cts = [(i-12)%22 for i in phases]
    N = 11
    radii = np.histogram(cts,N)[0]
    bottom = 0
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    width = (2*np.pi) / N
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(0.5*np.pi)
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(15)
    ax.set_xticks(np.pi/180. * np.linspace(0,  360, N, endpoint=False))
    ax.set_xticklabels(range(1,2*N+1,2))
    bars = ax.bar(theta, radii, width=width, bottom=bottom)
    # Use custom colors and opacity
    for r, bar in zip(radii, bars):
        bar.set_facecolor(plt.cm.jet(r / float(np.max(radii))))
        bar.set_alpha(0.8)
    plt.title(fname, y=1.08)
    plt.savefig(fname+'_phasehist.pdf')
    plt.close()


plt.scatter(prot['zscore_wt'].values,prot['zscore_csp'].values,c=np.abs(np.abs(prot['delta_p'].values)-12))
plt.scatter(rp['zscore_wt'].values,rp['zscore_rna'].values,c=np.abs(np.abs(rp['delta_p'].values)-12))




def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi

dpgmm = mixture.BayesianGaussianMixture(n_components=4, covariance_type='full').fit(X_train)

# display predicted scores by the model as a contour plot
x = np.linspace(-20., 30.)
y = np.linspace(-20., 40.)
X, Y = np.meshgrid(x, y)
XX = np.array([X.ravel(), Y.ravel()]).T
Z = dpgmm.predict_proba(XX)
Z = [Z[:,i].reshape(X.shape) for i in range(Z.shape[1])]

i=1
CS = plt.contour(X, Y, Z[i], locator=ticker.LogLocator(), levels=np.logspace(-3, 0, 10), cmap=colors[i])
CB = plt.colorbar(CS, shrink=0.8, extend='both')
plt.scatter(X_train[:, 0], X_train[:, 1], .8)

plt.title('Negative log-likelihood predicted by a GMM')
plt.axis('tight')
plt.ylim([-2,5])
plt.xlim([-2,5])
plt.show()




colors = ['Reds','Purples','Blues','Greys']

for i in range(len(Z)):
    CS = plt.contour(X, Y, Z[i], locator=ticker.LogLocator(), levels=np.logspace(-3, 0, 10), cmap=colors[i])

CB = plt.colorbar(CS, shrink=0.8, extend='both')
plt.scatter(X_train[:, 0], X_train[:, 1], .8)

plt.title('Negative log-likelihood predicted by a GMM')
plt.axis('tight')
plt.show()
