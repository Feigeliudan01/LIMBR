import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import numpy as np, pandas as pd, matplotlib.pyplot as plt; plt.rc('text', usetex=True)
import seaborn as sns; sns.set(style="white", context="talk")
from  matplotlib import colors, ticker, cm
from astroML.utils import convert_2D_cov
from astroML.plotting.tools import draw_ellipse
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=False)
from scipy import stats
from sklearn.cluster import KMeans
import scipy.cluster.hierarchy as sch
from sklearn.preprocessing import scale
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
circularlib=importr('circular')
circular=robjects.r('circular')
corcircular=robjects.r('cor.circular')
import pickle
from LIMBR import imputation

def plot_cormat(df,fname,ptitle):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    cors = np.corrcoef(df.values.T)
    z_min, z_max = -np.abs(cors).max(), np.abs(cors).max()
    pc = ax.pcolor(cors, cmap='seismic', vmin=z_min, vmax=z_max)
    for i in range(0,cors.shape[1],3):
        ax.axvline(x=i, linestyle='-', color='black', linewidth=0.7,    zorder=3)
    for i in range(0,cors.shape[0],3):
        ax.axhline(y=i, linestyle='-', color='black', linewidth=0.7,    zorder=3)
    tpts = [i*2 for i in range(1,(int(cors.shape[1]/3)+1))]
    ax.set_xticks([i+1.5 for i in range(0,cors.shape[1],3)])
    ax.set_xticklabels([(i-12)%22 for i in tpts])
    ax.set_yticks([i+1.5 for i in range(0,cors.shape[1],3)])
    ax.set_yticklabels([(i-12)%22 for i in tpts])
    ax.set_title(ptitle, fontsize=16)
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
    cb = fig.colorbar(pc, cax=cbar_ax)
    cb.set_label('Correlation', fontsize=12)
    plt.savefig(fname+'.pdf')
    plt.close()

def plot_cormat_prep(df,fname,ptitle):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    cors = np.corrcoef(df.values.T)
    z_min, z_max = -np.abs(cors).max(), np.abs(cors).max()
    pc = ax.pcolor(cors, cmap='seismic', vmin=z_min, vmax=z_max)
    for i in prep_breaks:
        ax.axvline(x=i, linestyle='-', color='black', linewidth=0.7,    zorder=3)
    for i in prep_breaks:
        ax.axhline(y=i, linestyle='-', color='black', linewidth=0.7,    zorder=3)
    ax.set_xticks([prep_breaks[i] + prep_diffs[i]/2 for i in range(len(prep_breaks))])
    ax.set_xticklabels([i+1 for i in range(16)])
    ax.set_yticks([prep_breaks[i] + prep_diffs[i]/2 for i in range(len(prep_breaks))])
    ax.set_yticklabels([i+1 for i in range(16)])
    ax.set_title(ptitle, fontsize=16)
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
    cb = fig.colorbar(pc, cax=cbar_ax)
    cb.set_label('Correlation', fontsize=12)
    plt.savefig(fname+'.pdf')
    plt.close()

def scale_df(df):
    sc = df.apply(lambda x: x-np.mean(x),axis=1)
    sc = sc.apply(lambda x: x/np.std(x),axis=1)
    return sc

def reorder_heatmap(arr):
    Y = sch.linkage(arr, method='complete', metric='euclidean')
    Z = sch.dendrogram(Y, orientation='left',no_plot=True)
    idx = Z['leaves']
    return arr[idx,:]

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
    plt.close()


def autocorr(l,shift):
    return np.dot(l, np.roll(l, shift)) / np.dot(l, l)

def autocorr(x,shift):
    result = np.correlate(x, x, mode='full')
    return result[result.size//2:][shift]

def get_tpoints(df):
    tpoints = [i.replace('CT','') for i in df.columns.values]
    tpoints = [int(i.split('.')[0]) if '.' in i else int(i.split('_')[0]) for i in tpoints]
    return np.asarray(tpoints)

def circ_cor(df):
    per = 6
    cors = []
    tpoints = get_tpoints(df)
    for row in df.values:
        ave = []
        #might eventually need to account for case where all replicates of a timepoint are missing (in this case the experiment is probably irreparably broken anyway though)
        for k in set(tpoints):
            ave.append((np.mean([row[i] for i, j in enumerate(tpoints) if j == k])))
        cors.append(autocorr(ave,(per*2))-autocorr(ave,per))
    return np.asarray(cors)

def get_fstats(ts,cs,ctypes):
    results = {}
    results2 = {}
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
        tempicc = {}
        for k, v in fortest.items():
            tempdict[k] = stats.f_oneway(*v)[1]
            tempicc[k] = icc(v)
        results[class_val] = tempdict
        results2[class_val] = tempicc
    return pd.DataFrame.from_dict(results), pd.DataFrame.from_dict(results2)


def gen_biasmap(df1,df2,fname,ptitle):
    df1.sort_index(ascending=False,inplace=True)
    df2.sort_index(ascending=False,inplace=True)
    fig, ax = plt.subplots(nrows=1, ncols=1)
    pc = ax.pcolor(df1.values, cmap='viridis')
    ax.set_xticks([i+.5 for i in range(len(df1.columns.values))])
    ax.set_xticklabels(df1.columns.values)
    ax.set_yticks([i+.5 for i in range(len(df1))])
    ylabs = [i for i in range(1,len(df1)+1)]
    ylabs.reverse()
    ax.set_yticklabels(ylabs)
    ax.xaxis.tick_top()
    ax.set_ylabel('Bias Trend')
    ax.set_title(ptitle, fontsize=16, y=1.08)
    for y in range(df2.values.shape[0]):
        for x in range(df2.values.shape[1]):
            if df2.values[y,x] < .001/df2.values.shape[1]:
                plt.text(x + 0.5, y + 0.5, '*', horizontalalignment='center',verticalalignment='center', fontsize=16)
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
    cb = fig.colorbar(pc, cax=cbar_ax)
    cb.set_label('ICC', fontsize=12)
    plt.savefig(fname+'.pdf')
    plt.close()

def icc(arr):
    m = np.mean([np.mean(i) for i in arr])
    alphas = [np.mean(i) for i in arr]-m
    es = [list(arr[i]-m-alphas[i]) for i in range(len(arr))]
    es = [item for sublist in es for item in sublist]
    return (np.std(alphas)**2)/(np.std(alphas)**2+np.std(es)**2)

def scale_peps(df):
    sc = df.apply(lambda x: x-np.mean(x),axis=0)
    sc = sc.apply(lambda x: x/np.std(x),axis=0)
    return sc

def plot_cormat_pep(arr,fname,ptitle):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    #z_min, z_max = -np.abs(arr).max(), np.abs(arr).max()
    z_min, z_max = arr.min(), arr.max()
    pc = ax.pcolor(arr, cmap='seismic', vmin=z_min, vmax=z_max)
    ax.set_xticks([i+.5 for i in range(0,arr.shape[1],1)])
    ax.set_yticks([i+.5 for i in range(0,arr.shape[1],1)])
    ax.set_xticklabels(range(1,arr.shape[1]+1))
    ax.set_yticklabels(range(1,arr.shape[1]+1))
    ax.set_title(ptitle, fontsize=16)
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
    cb = fig.colorbar(pc, cax=cbar_ax)
    cb.set_label('Correlation', fontsize=12)
    plt.savefig(fname+'.pdf')
    plt.close()

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
    ax.text(0.85, 0.05, r'$\rho_{cc}$ = %f'"\n"r'p = %e'%(float(circorr[0].r_repr()), float(circorr[2].r_repr())), transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=dict(facecolor='white'))
    plt.tight_layout()
    plt.savefig('output/figs/WT_CSP_phasehist.pdf')
    plt.close()

wt_p = pd.read_csv('output/actual/wt_lowess_normed__jtkout_GammaP.txt',sep='\t',index_col=0)
csp_p = pd.read_csv('output/actual/csp_lowess_normed__jtkout_GammaP.txt',sep='\t',index_col=0)
rna_p = pd.read_csv('output/actual/rna_lowess_normed__jtkout_GammaP.txt',sep='\t',index_col=0)
wt_old_p = pd.read_csv('output/actual/wt_old_normed__jtkout_GammaP.txt',sep='\t',index_col=0)
csp_old_p = pd.read_csv('output/actual/csp_old_normed__jtkout_GammaP.txt',sep='\t',index_col=0)
rna_old_p = pd.read_csv('output/actual/rna_old_normed__jtkout_GammaP.txt',sep='\t',index_col=0)

wt = pd.read_csv('output/actual/wt_lowess_normed.txt',sep='\t',index_col=0)
csp = pd.read_csv('output/actual/csp_lowess_normed.txt',sep='\t',index_col=0)
rna = pd.read_csv('output/actual/rna_lowess_normed.txt',sep='\t',index_col=0)
wt_old = pd.read_csv('output/actual/wt_old_normed.txt',sep='\t',index_col=0)
csp_old = pd.read_csv('output/actual/csp_old_normed.txt',sep='\t',index_col=0)
rna_old = pd.read_csv('output/actual/rna_old_normed.txt',sep='\t',index_col=0)


#Figure 1A
wt_class = pd.read_csv('data/wt_classes.csv')
csp_class = pd.read_csv('data/csp_classes.csv')

plot_cormat(scale_df(wt_old),'output/figs/initial_cormat','Correlation Matrix of Traditionally Normalized WT Samples Ordered by Timepoint')
new_cols = [i.replace('_1','') for i in wt_class[~wt_class['Sample ID'].isin(['CT00_3', 'CT00_2', 'CT00_1'])]['Sample ID'].values]
new_cols = [i.replace('_2','.1') for i in new_cols]
new_cols = [i.replace('_3','.2') for i in new_cols]
#new_cols = [i for i in new_cols if i not in ['CT00.2', 'CT00.1', 'CT00']]
prep = wt_class[~wt_class['Sample ID'].isin(['CT00_3', 'CT00_2', 'CT00_1'])]['Prep Set'].values
prep_breaks = []
curr_class = prep[0]
for i in range(len(prep)):
    if prep[i] != curr_class:
        prep_breaks.append(i)
    curr_class = prep[i]

prep_diffs = [4] + [prep_breaks[i+1]-prep_breaks[i] for i in range(len(prep_breaks)-1)] +[2]
prep_breaks = [0] + prep_breaks

plot_cormat_prep(scale_df(wt_old[new_cols]),'output/figs/initial_cormat_prep','Correlation Matrix of Traditionally Normalized WT Samples Ordered by Prep Set')


#Figure 1B
rawdata = pd.read_csv('data/Formatted_input.txt',sep='\t')
dist = rawdata.isnull().sum(axis=1).values.tolist()
density = stats.gaussian_kde(dist)
xs = np.linspace(0,max(dist),200)
density.covariance_factor = lambda : .1
density._compute_covariance()
values, base = np.histogram(dist, bins=100)
cumulative = np.cumsum(values)

fig, ax1 = plt.subplots()
ax1.plot(xs,density(xs),'b')
ax1.set_xlabel('Number of Missing Values')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('Density of Peptides', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(base[:-1], cumulative, c='r')
ax2.set_ylabel('Cumulative Number of Peptides \n at or Below Missingness Level', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.savefig('output/figs/missingness.pdf')
plt.close()

#Figure 2A
simres = pd.read_csv('output/simdata/simdata_mb.csv',index_col=0)
plt.scatter(simres['Base_FPR'].values,simres['Base_TPR'].values,marker='x', label='Baseline classification without noise',color='g')
plt.scatter(simres['Circ_FPR'].values,simres['Circ_TPR'].values,marker='x', label='Circadian denoised classification',color='b')
plt.scatter(simres['Noise_FPR'].values,simres['Noise_TPR'].values,marker='x', label='Classification with noise',color='r')
plt.axvline(x=0.05,color='black',ls='dashed')
plt.legend(loc="lower right")
plt.xlabel('1 - Specificity')
plt.ylabel('Sensitivity')
plt.xlim(0,1)
plt.ylim(0,1)
plt.savefig('output/figs/sim_ROC_mb.pdf')
plt.close()

#Figure 2B
simres0 = pd.read_csv('output/simdata/simdata_mb.csv')
melted_simres = pd.melt(simres0, id_vars=['Iteration'], value_vars=[r'Base_auc', r'Noise_auc',r'Circ_auc'])
ax = sns.boxplot(x='variable',y='value',data=melted_simres)
ax = sns.swarmplot(x='variable',y='value',data=melted_simres)
plt.savefig('output/figs/sim_ROC_auc_mb.pdf')
plt.close()

#Figure 2C
key = pd.read_csv('output/simdata/mb_simulated_data_key_1.txt',sep='\t',index_col=0)
simcirc_base = pd.read_csv('output/simdata/mb_simulated_data_baseline_1.txt',sep='\t',index_col=0)[key['circ']==1]
simcirc_noise = pd.read_csv('output/simdata/mb_simulated_data_with_noise_1.txt',sep='\t',index_col=0)[key['circ']==1]
simcirc_denoise = pd.read_csv('output/simdata/mb_denoised_circ_lowess_1.txt',sep='\t',index_col=0)[key['circ']==1]

plot_heatmap(simcirc_base.values,r'output/figs/simulated_circ_genes_mb','Expression as a Function of Timepoint in Simulated Circadian Genes')
plot_heatmap(simcirc_noise.values,r'output/figs/simulated_noise_circ_genes_mb','Expression as a Function of Timepoint in Simulated Circadian Genes with Noise')
plot_heatmap(simcirc_denoise.values,r'output/figs/simulated_denoise_circ_genes_mb','Expression as a Function of Timepoint in Denoised Simulated Circadian Genes')

#Figure S1
acorr = key.copy()
acorr['cors'] = circ_cor(pd.read_csv('output/simdata/mb_simulated_data_baseline_1.txt',sep='\t',index_col=0))
acorr['Circadian'] = acorr['circ'].apply(lambda x: 'Circ' if x == 1 else 'Non-Circ')
acorr['included'] = [(i<(np.percentile(acorr['cors'].values,25))) for i in acorr['cors'].values]
plt.axhline(y=np.percentile(acorr['cors'].values,25),color='black',ls='dashed')
sns.swarmplot(x='circ', y='cors', hue='included', data = acorr)
plt.ylabel(r'$\Delta$ Autocorrelation')
plt.savefig('output/figs/acorr_val_mb.pdf')
plt.close()

#Fig 3 B and C
def get_point(fname_jtk, fname_key):
    jtk = pd.read_csv(fname_jtk,sep='\t', index_col=0)
    key = pd.read_csv(fname_key,sep='\t', index_col=0)
    merged = pd.merge(jtk, key,left_index=True,right_index=True)
    return ((1-len(merged[(merged['GammaP']>.05)&(merged['circ']==0)])/len(merged[merged['circ']==0]))),(len(merged[(merged['GammaP']<.05)&(merged['circ']==1)])/len(merged[merged['circ']==1]))

def get_auc(fname_jtk, fname_key):
    jtk = pd.read_csv(fname_jtk,sep='\t', index_col=0)
    key = pd.read_csv(fname_key,sep='\t', index_col=0)
    merged = pd.merge(jtk, key,left_index=True,right_index=True)
    scoresn = np.asarray([1-i for i in merged['GammaP'].values.astype(float)])
    yn = merged['circ'].values.astype(int)
    fprn, tprn, thresholdsn = metrics.roc_curve(yn, scoresn)
    return metrics.auc(fprn, tprn)

baseline = []
noise = []
eigenms = []
limbr = []

for i in range(1,20):
    baseline.append(get_point('output/simdata/mb_simulated_data_baseline_'+str(i)+'__jtkout_GammaP.txt','output/simdata/mb_simulated_data_key_'+str(i)+'.txt'))
    noise.append(get_point('output/simdata/simulated_data_with_noise_'+str(i)+'__jtkout_GammaP.txt','output/simdata/mb_simulated_data_key_'+str(i)+'.txt'))
    eigenms.append(get_point('output/simdata/mb_simdata_eigenMS_'+str(i)+'__jtkout_GammaP.txt','output/simdata/mb_simulated_data_key_'+str(i)+'.txt'))
    limbr.append(get_point('output/simdata/mb_denoised_circ_lowess_'+str(i)+'__jtkout_GammaP.txt','output/simdata/mb_simulated_data_key_'+str(i)+'.txt'))

plt.scatter(*zip(*baseline),color='b',marker='x', label='Baseline ROC curve')

plt.scatter(*zip(*noise),color='r',marker='x', label='Noise ROC curve')

plt.scatter(*zip(*eigenms),color='y',marker='x', label='EigenMS ROC curve')


plt.scatter(*zip(*limbr),color='g',marker='x', label='LIMBR ROC curve')


plt.axvline(x=0.05,color='black',ls='dashed')
plt.legend(loc="lower right")
plt.axis([0, 1, 0, 1])
plt.xlabel('1 - Specificity')
plt.ylabel('Sensitivity')
plt.savefig('output/figs/sim_ROC_mb.pdf')
plt.close()

baseline_auc = []
noise_auc = []
eigenms_auc = []
limbr_auc = []

for i in range(1,20):
    baseline_auc.append(get_auc('output/simdata/mb_simulated_data_baseline_'+str(i)+'__jtkout_GammaP.txt','output/simdata/mb_simulated_data_key_'+str(i)+'.txt'))
    noise_auc.append(get_auc('output/simdata/simulated_data_with_noise_'+str(i)+'__jtkout_GammaP.txt','output/simdata/mb_simulated_data_key_'+str(i)+'.txt'))
    eigenms_auc.append(get_auc('output/simdata/mb_simdata_eigenMS_'+str(i)+'__jtkout_GammaP.txt','output/simdata/mb_simulated_data_key_'+str(i)+'.txt'))
    limbr_auc.append(get_auc('output/simdata/mb_denoised_circ_lowess_'+str(i)+'__jtkout_GammaP.txt','output/simdata/mb_simulated_data_key_'+str(i)+'.txt'))

auc = pd.DataFrame({'baseline' : baseline_auc, 'noise' : noise_auc, 'eigenms' : eigenms_auc, 'LIMBR' : limbr_auc})

auc = pd.melt(auc)

auc.columns = ['Processing','Area Under the Curve']

sns.boxplot(x="Processing", y="Area Under the Curve", data=auc, order=['baseline','noise','eigenms','LIMBR'], palette=['b','r','y','g'])

plt.savefig('output/figs/sim_AUC_mb.pdf')
plt.close()

#Figure 4A

wt_trends = pd.read_csv('output/actual/wt_lowess_normed_trends.txt',sep='\t')
wt_trends = wt_trends.drop('Unnamed: 0',axis=1)
wt_trends['Trend'] = range(1,len(wt_trends)+1)
wt_trends = pd.melt(wt_trends, id_vars=['Trend'], value_vars=wt_trends.columns.values.tolist()[:-1])


wt_bias_f, wt_bias_icc = get_fstats(wt_trends,wt_class,['TMT Set', 'Prep Set'])

gen_biasmap(wt_bias_icc,wt_bias_f,'output/figs/wt_bcor',r'Correlation of Bias Trends with Groupings in Experimental Procedure (WT)')

#Figure 4B
csp_trends = pd.read_csv('output/actual/csp_lowess_normed_trends.txt',sep='\t')
csp_trends = csp_trends.drop('Unnamed: 0',axis=1)
csp_trends['Trend'] = range(1,len(csp_trends)+1)
csp_trends = pd.melt(csp_trends, id_vars=['Trend'], value_vars=csp_trends.columns.values.tolist()[:-1])

csp_bias_f, csp_bias_icc = get_fstats(csp_trends,csp_class,['TMT Set', 'Prep Set'])

gen_biasmap(csp_bias_icc,csp_bias_f,'output/figs/csp_bcor',r'Correlation of Bias Trends with Groupings in Experimental Procedure ($\Delta$CSP-1)')

#Figure 4C
wt_pepts =  pd.read_csv('output/actual/wt_lowess_normed_pep_bias.txt',sep='\t')
csp_pepts =  pd.read_csv('output/actual/csp_lowess_normed_pep_bias.txt',sep='\t')
merge_pepts = pd.merge(wt_pepts.drop('Protein',axis=1), csp_pepts.drop('Protein',axis=1), on='Peptide').set_index('Peptide')
merge_pepts = scale(merge_pepts)

pep_corr = np.corrcoef(merge_pepts.T)[(len(wt_pepts.columns)-2):((len(wt_pepts.columns)-2)+5),0:5]

plot_cormat_pep(np.abs(pep_corr),'output/figs/peptide_btrend_cormat','Correlation Matrix of Peptide Level Bias Trends')

#Figure 5A
plot_cormat(scale_df(wt),'output/figs/final_cormat','Correlation Matrix of LIMBR Normalized WT Samples Ordered by Timepoint')

prot = pd.merge(wt_p[['GammaP','Phase']], csp_p[['GammaP','Phase']], right_index=True, left_index=True,suffixes=('_wt','_csp'))
prot['zscore_wt'] = stats.norm.ppf(1-prot['GammaP_wt'])
prot['zscore_csp'] = stats.norm.ppf(1-prot['GammaP_csp'])
prot['delta_z'] = prot['zscore_wt'] - prot['zscore_csp']
prot['delta_p'] = prot['Phase_csp'] - prot['Phase_wt']

plot_cormat(scale_df(wt[wt.index.isin(prot[prot['GammaP_wt']<.05].index.tolist())]),'output/figs/final_cormat_circ','Correlation Matrix of LIMBR Normalized WT Circadian Proteins Ordered by Timepoint')

#Figure 5B

defcirc = prot[(prot['GammaP_wt']<.05)&(prot['GammaP_csp']<.05)]

cpw = circular(robjects.FloatVector([j*2*np.pi/22 for j in defcirc['Phase_wt'].values]),units="radians", zero=0, rotation='clock')
cpc = circular(robjects.FloatVector([j*2*np.pi/22 for j in defcirc['Phase_csp'].values]),units="radians", zero=0, rotation='clock')

circorr = corcircular(cpw,cpc,test='TRUE')

phase_hist(defcirc['delta_p'].values,r'Histogram of Phase Lags Between'"\n"r'Expression in WT and $\Delta$CSP-1 Genotypes'"\n"r'For Proteins Circadian in Both Genotypes')


#Figure 1d
to_impute = imputation.imputable("data/Formatted_input.txt",.3)
to_impute.deduplicate()
miss = []
blocks = []
for i in range(1,20):
    miss.append(len((to_impute.data[to_impute.data.columns[0:(i*10)]]).dropna()))
    blocks.append(i)

fig, ax1 = plt.subplots()
ax1.plot(blocks,miss, color='b')[0]
ax1.set_ylabel('Number of peptides without any missing values', color='k')
ax1.set_xlabel('Number of MS runs', color='k')
ax1.set_xticks(range(1,20))
plt.savefig('output/figs/missingness_vs_experiments.pdf')
plt.close()

#Figure ???
#ax = plt.subplot(111)
#hb = ax.hexbin(prot['zscore_wt'].values, prot['zscore_csp'].values, gridsize=48, extent=(-2,4,-2,4))
#ax.set_xlabel(r"WT Z Score", fontsize=12)
#ax.set_ylabel(r'$\Delta$CSP-1 Z Score', fontsize=12)
#cb = fig.colorbar(hb, ax=ax, ticks=np.linspace(0,20,11))
#cb.set_label('count', fontsize=12)
#ax.set_title('Distribution of Z Scores from eJTK of Proteomics Data by Genotype', fontsize=16)
#ax.text(0.015, 0.985, "Pearsons r = %f, p = %e"%(r, p_value), transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=dict(facecolor='white'))
#plt.savefig('prot_hex.pdf')
#plt.close()

#Figure ???
#sns.jointplot('zscore_wt','zscore_csp',data=prot,kind='kde',joint_kws=dict(bw=.1),marginal_kws=dict(bw=.1))

#Figure ???
