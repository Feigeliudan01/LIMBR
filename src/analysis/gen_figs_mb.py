import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import numpy as np, pandas as pd, matplotlib.pyplot as plt; plt.rc('text', usetex=True)
import seaborn as sns; sns.set(style="white", context="talk")
from  matplotlib import colors, ticker, cm
from astroML.utils import convert_2D_cov
from astroML.plotting.tools import draw_ellipse
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=True)
from scipy import stats
from sklearn.cluster import KMeans
import scipy.cluster.hierarchy as sch
from sklearn.preprocessing import scale
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
circularlib=importr('circular')
circular=robjects.r('circular')
corcircular=robjects.r('cor.circular')

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
    tpoints = [int(i.split('.')[0]) for i in tpoints]
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
    z_min, z_max = -np.abs(arr).max(), np.abs(arr).max()
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
acorr['cors'] = circ_cor(pd.read_csv('output/simdata/_mbsimulated_data_baseline_1.txt',sep='\t',index_col=0))
acorr['Circadian'] = acorr['circ'].apply(lambda x: 'Circ' if x == 1 else 'Non-Circ')
acorr['included'] = [(i<(np.percentile(acorr['cors'].values,25))) for i in acorr['cors'].values]
plt.axhline(y=np.percentile(acorr['cors'].values,25),color='black',ls='dashed')
sns.swarmplot(x='circ', y='cors', hue='included', data = acorr)
plt.ylabel(r'$\Delta$ Autocorrelation')
plt.savefig('output/figs/acorr_val_mb.pdf')
plt.close()
