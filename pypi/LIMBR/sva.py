import numpy as np
import pandas as pd
import os
import time
import scipy.stats as stats
from numpy.linalg import svd, lstsq
from sklearn.decomposition import PCA
from scipy.stats import linregress, f_oneway
import itertools
import sys
from statsmodels.nonparametric.smoothers_lowess import lowess
from tqdm import tqdm
from sklearn.preprocessing import scale
from sklearn.neighbors import NearestNeighbors
import math
import json
from ctypes import c_int
import pickle
from multiprocess import Pool, current_process, Manager
from functools import partial
from sklearn import preprocessing

class sva:

    def __init__(self, filename,design,data_type,blocks=None,pool=None):
        """
        imports data and initializes an sva object

        takes a file from one of two data types protein ('p') which has two index columns or rna ('r') which has only one.

        opens a pickled file matching pooled controls to corresponding samples
        """
        np.random.seed(4574)
        self.data_type = str(data_type)
        if self.data_type == 'p':
            self.raw_data = pd.read_csv(filename,sep='\t').set_index(['Peptide','Protein'])
        if self.data_type == 'r':
            self.raw_data = pd.read_csv(filename,sep='\t').set_index('#')
        self.designtype = str(design)
        if self.designtype == 'b':
            self.block_design = pickle.load( open( blocks, "rb" ) )
        if pool != None:
            self.norm_map = pickle.load( open( pool, "rb" ) )
        elif pool == None:
            self.norm_map = None
        self.notdone = True

    def pool_normalize(self):
        """
        performs pool normalization on an sva object using the raw_data and norm_map if pooled controls were used, then performs quantile normalization.

        if no pooled controls were used as for RNAseq experiments, only quantile normalization is performed.

	after normalization, each row is scaled.
        """

        def pool_normalize(df,dmap):
            """
            divides peptide abundances of each sample by corresponding pooled control abundance

	    takes a dataframe from raw_data and a map of pooled controls from norm_map

            returns a pool normalized dataframe
            """
            newdf = pd.DataFrame(index=df.index)
            for column in df.columns.values:
                newdf[column] = df[column].div(df['pool_'+'%02d' % dmap[column]],axis='index')
            nonpool = [i for i in newdf.columns if 'pool' not in i]
            newdf = newdf[nonpool]
            return newdf

        def qnorm(df):
            """
            quantile normalizes data by columns

	    takes a dataframe either raw_data or data_pnorm

	    returns a quantile normalized dataframe
            """
            ref = pd.concat([df[col].sort_values().reset_index(drop=True) for col in df], axis=1, ignore_index=True).mean(axis=1).values
            for i in range(0,len(df.columns)):
                df = df.sort_values(df.columns[i])
                df[df.columns[i]] = ref
            return df.sort_index()

        if (self.data_type == 'r') or (self.norm_map == None):
            self.data = qnorm(self.raw_data)
            self.scaler = preprocessing.StandardScaler().fit(self.data.values.T)
            self.data = pd.DataFrame(self.scaler.transform(self.data.values.T).T,columns=self.data.columns,index=self.data.index)
        else:
            self.data_pnorm = pool_normalize(self.raw_data,self.norm_map)
            self.data_pnorm = self.data_pnorm.replace([np.inf, -np.inf], np.nan)
            self.data_pnorm = self.data_pnorm.dropna()
            self.data_pnorm = self.data_pnorm.sort_index(axis=1)
            self.data_pnorm = qnorm(self.data_pnorm)
            self.scaler = preprocessing.StandardScaler().fit(self.data_pnorm.values.T)
            self.data = pd.DataFrame(self.scaler.transform(self.data_pnorm.values.T).T,columns=self.data_pnorm.columns,index=self.data_pnorm.index)


    def get_tpoints(self):
        """
        extracts timepoints from header of data
        """
        tpoints = [i.replace('CT','') for i in self.data.columns.values]
        tpoints = [int(i.split('_')[0]) for i in tpoints]
        #deprecated splitting for alternative header syntax
        #tpoints = [int(i.split('.')[0]) for i in tpoints]
        self.tpoints = np.asarray(tpoints)

    def prim_cor(self):
        """
        calculates correlation for each row against the primary variable of interest.
        primary variable of interest is defined based on the designtype.

        for circadian ('c') this is the difference between the autocorrelation at the expected period and one half that period
        for general timecourse ('l') this is the goodness of fit of a lowess model
        for a block design ('b') this is the ANOVA f statistic
        """
        def circ_cor():
            """
            calculates rough estimate of circadianness based on autocorrelation differences

            given a period fixed here at 12 samples the autocorrelation is calculated at 12 and 6 samples shift and the difference indicates the degree of circadian signal
            """
            def autocorr(l,shift):
                """calculates autocorrelation of a series l at a given shift"""
                return np.dot(l, np.roll(l, shift)) / np.dot(l, l)

            per = 12
            cors = []
            for row in tqdm(self.data.values):
                ave = []
                #might eventually need to account for case where all replicates of a timepoint are missing (in this case the experiment is probably irreparably broken anyway though)
                for k in set(self.tpoints):
                    ave.append((np.mean([row[i] for i, j in enumerate(self.tpoints) if j == k])*1000000))
                cors.append((autocorr(ave,per) - autocorr(ave,(per//2))))
            self.cors = np.asarray(cors)

        def l_cor():
            """
            calculates how well a row of data fits expectations for a timecourse

            this is defined as the goodness of fit of a lowess model for the row
            """
            cors = []
            for row in tqdm(self.data.values):
                ys = lowess(row, self.tpoints, it=1)[:,1]
                cors.append(-sum((row - ys)**2))
            self.cors = np.asarray(cors)

        def block_cor():
            """
            calculates how well a row of data fits expectations for a block study design

            this is defined as the ANOVA f statistic which compares between and within block variability
            """
            cors = []
            for row in tqdm(self.data.values):
                blist = []
                for k in set(self.block_design):
                    blist.append(([row[i] for i, j in enumerate(self.block_design) if j == k]))
                cors.append(f_oneway(*blist)[0])
            self.cors = np.asarray(cors)

        if self.designtype == 'c':
            circ_cor()
        elif self.designtype == 'b':
            block_cor()
        elif self.designtype == 'l':
            l_cor()


    def reduce(self,percsub):
        """
        reduces the data based on the correlation to primary variable of interest calculated in prim_cor

        takes cors and a percentage cutoff and removes that percentage of data most correlated with variable of interest
        """
        percsub = float(percsub)
        uncor = [(i<(np.percentile(self.cors,percsub))) for i in self.cors]
        self.data_reduced = self.data[uncor]


    def get_res(self,in_arr):
        """
        calculates model residuals from which to learn batch effects

        this is done either with a lowess model for timecourse designs or with the block mean for block designs
        """
        def get_l_res(arr):
            """calculates residuals from a lowess model for timecourse designs"""
            res = []
            for row in arr:
                ys = lowess(row, self.tpoints, it=1)[:,1]
                res.append(row - ys)
            return np.array(res)

        def get_b_res(arr):
            """calculates residuals from the block mean for block based designs"""
            m = {}
            for v in set(self.block_design):
                indices = [i for i, x in enumerate(self.block_design) if x == v]
                m[v] = np.mean(arr[:,indices],axis=1)
            ma = np.zeros(np.shape(arr))
            for i in tqdm(range(len(self.block_design)), desc='get block residuals 2', leave=False):
                ma[:,i]=m[self.block_design[i]]
            return np.subtract(arr,ma)
        if self.designtype == 'c':
            return get_l_res(in_arr)
        elif self.designtype == 'b':
            return get_b_res(in_arr)
        elif self.designtype == 'l':
            return get_l_res(in_arr)

    #defined seperately for reuse in perm_test
    def set_res(self):
        """calls get_res()"""
        self.res = self.get_res(self.data_reduced.values)

    def get_tks(self,arr):
        """calculates the fraction of variance for each row of the residual matrix explained by each principle component with PCA"""
        pca = PCA(svd_solver='randomized',random_state=4574)
        pca.fit(arr)
        return pca.explained_variance_ratio_

    #defined seperately for reuse in perm_test
    def set_tks(self):
        """calls get_tks()"""
        self.tks = self.get_tks(self.res)

    def perm_test(self,nperm,npr):
        """
        performs permutation testing on residual matrix PCA

        multiprocessing with number of processors npr

        runs single_it number of times specified by nperm and returns how many times tkstar was bigger than tk for all tks
        """
        mgr = Manager()
        output = mgr.list()
        def single_it(rseed):
            """single iteration of permutation testing, permutes residual matrix, calculates new tks for permuted matrix and compares to original tks"""
            rstate = np.random.RandomState(rseed*100)
            rstar = np.copy(self.res)
            out = np.zeros(len(self.tks))
            for i in range(rstar.shape[0]):
                rstate.shuffle(rstar[i,:])
            resstar = self.get_res(rstar)
            tkstar = self.get_tks(resstar)
            for m in range(len(self.tks)):
                if tkstar[m] > self.tks[m]:
                    out[m] += 1
            return out
        l = mgr.Lock()
        with Pool(int(npr)) as pool:
            pbar = tqdm(total=int(nperm), desc='permuting', position=0, smoothing=0)
            imap_it = pool.imap_unordered(single_it, range(int(nperm)))
            for x in imap_it:
                pbar.update(1)
                with l:
                    output.append(x)
        pbar.close()
        pool.close()
        pool.join()
        self.sigs = np.sum(np.asarray(output), axis=0)/float(nperm)
        time.sleep(40)

    def eig_reg(self,alpha):
        """regresses eigentrends (batch effects) against the reduced dataset calculating p values for each row being associated with that trend"""
        alpha = float(alpha)
        U, s, V = np.linalg.svd(self.res)
        sig = V.T[:,:len([i for i in itertools.takewhile(lambda x: x < alpha, self.sigs.copy())])]
        pvals = []
        if len(sig)>0:
            for trend in tqdm(sig.T.copy()):
                temp = []
                for row in self.data_reduced.values.copy():
                    slope, intercept, r_value, p_value, std_err = linregress(row,trend)
                    temp.append(p_value)
                pvals.append(temp)
            self.ps =  pvals
        else:
            print('No Significant Trends')

    def subset_svd(self,lam):
        """
        performs SVD on the subset of rows associated with each batch effect

        after subsetting the reduced data, SVD is performed on this matrix and the true batch effect is taken as the right singular vector most correlated with the initial batch effect from np.linalg.svd(self.res)
        """

        def est_pi_naught(probs_naught,lam):
            """estimates background distribution of p values"""
            return len([i for i in probs_naught if i > lam])/(len(probs_naught)*(1-lam))

        def est_pi_sig(probs_sig,l):
            """finds p value pi_sig as cutoff for rows associated with batch effect"""
            pi_0 = est_pi_naught(probs_sig,l)
            if pi_0 > 1:
                return 'nan'
            sp = np.sort(probs_sig)
            return sp[int(np.floor((1-pi_0)*len(probs_sig))-1)]

        pt, _, bt = np.linalg.svd(self.res)
        trends = []
        pep_trends = []
        for j, entry in enumerate(tqdm(self.ps)):
            sub = []
            thresh = est_pi_sig(entry,lam)
            if thresh == 'nan':
                self.ts = trends
                self.pepts = pep_trends
                return
            for i in range(len(entry)):
                if entry[i] < thresh:
                    sub.append(self.data_reduced.values[i])
            U, s, V = np.linalg.svd(sub)
            temp = []
            for trend in V:
                _, _, _, p_value, _ = linregress(bt[j],trend)
                temp.append(p_value)
            trends.append(V.T[:,np.argmin(temp)])
            pep_trends.append(pt[:,j])
        self.pepts = pep_trends
        self.ts = trends


    def normalize(self,outname):
        """
        creates diagnostic files, normalizes data based on calculated batch effects and outputs final processed dataset
        """
        pd.DataFrame(self.ts,columns=self.data.columns).to_csv(outname.split('.txt')[0]+'_trends.txt',sep='\t')
        pd.DataFrame(self.sigs).to_csv(outname.split('.txt')[0]+'_perms.txt',sep='\t')
        pd.DataFrame(self.tks).to_csv(outname.split('.txt')[0]+'_tks.txt',sep='\t')
        pd.DataFrame(np.asarray(self.pepts).T,index=self.data_reduced.index).to_csv(outname.split('.txt')[0]+'_pep_bias.txt',sep='\t')
        fin_res = np.dot(np.linalg.lstsq(np.asarray(self.ts).T,self.data.values.T)[0].T,np.asarray(self.ts))
        self.svd_norm = self.scaler.inverse_transform((self.data.values - fin_res).T).T
        self.svd_norm = pd.DataFrame(self.svd_norm,index=self.data.index,columns=self.data.columns)
        if self.data_type == 'p':
            self.svd_norm = self.svd_norm.groupby(level='Protein').mean()
        self.svd_norm.index.names = ['#']
        self.svd_norm.to_csv(outname,sep='\t')