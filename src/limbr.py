from random import *
import numpy as np
import pandas as pd
import os
import scipy.stats as stats
from numpy.linalg import svd, lstsq
import timeit
from sklearn.decomposition import RandomizedPCA
from scipy.stats import linregress
import itertools
import sys
import getopt
from pylab import dot, floor
from statsmodels.nonparametric.smoothers_lowess import lowess
from tqdm import tqdm
from sklearn.preprocessing import scale
from sklearn.neighbors import NearestNeighbors
import multiprocessing, logging
import math
import json
import time
from ctypes import c_int

class imputable:

    def __init__(self, filename, missingness):
        self.data = pd.read_csv(filename,sep='\t')
        self.miss = float(missingness)
        self.pats = {}
        self.notdone = True

    def deduplicate(self):
        if self.data[self.data.columns.values[1]][0][-2] == "T":
            self.data[self.data.columns.values[1]] = self.data[self.data.columns.values[1]].apply(lambda x: x.split('T')[0])
            self.data = self.data.groupby(['Peptide','Protein']).mean()
        todrop = []
        for name, group in tqdm(self.data.groupby(level='Peptide')):
            if len(group) > 1:
                todrop.append(name)
        self.data = self.data.drop(todrop)

    def drop_missing(self):
        self.miss = np.rint(len(self.data.columns)*self.miss)
        self.data = self.data[self.data.isnull().sum(axis=1)<=self.miss]

    def impute(self,outname):

        def match_pat(l,i):
            #make vector of missing/non-missing
            l = "".join(np.isnan(l).astype(int).astype(str))
            if l not in self.pats.keys():
                self.pats[l] = [i]
            else:
                self.pats[l].append(i)

        def get_patterns(arr):
            for ind, val in enumerate(arr):
                match_pat(val,ind)

        def sub_imputer(inds,pattern,origarr,comparr):
            #drop missing columns given missingness pattern
            newarr = comparr[:,~np.array(list(pattern)).astype(bool)]
            #fit nearest neighbors
            nbrs = NearestNeighbors(n_neighbors=10).fit(newarr)
            outa = []
            #iterate over rows matching missingness pattern
            for rowind, row in enumerate(origarr[inds]):
                outl = []
                #get indexes of given rows nearest neighbors
                indexes = nbrs.kneighbors([origarr[inds[rowind],~np.array(list(pattern)).astype(bool)]],return_distance=False)
                #get array of nearest neighbors
                means = np.mean(comparr[indexes[0][1:]], axis=0)
                #iterate over entries in each row
                for ind, v in enumerate(row):
                    if not np.isnan(v):
                        outl.append(v)
                    else:
                        outl.append(means[ind])
                outa.append(outl)
            return outa

        def imputer(origarr, comparr):
            outdict = {}
            for k in tqdm(self.pats.keys()):
                temparr = sub_imputer(self.pats[k],k, origarr,comparr)
                for ind, v in enumerate(temparr):
                    outdict[self.pats[k][ind]] = v
            return outdict


        datavals = self.data.values
        comparr = datavals[~np.isnan(datavals).any(axis=1)]
        get_patterns(datavals)
        out = imputer(datavals, comparr)
        meld = pd.DataFrame.from_dict(out,orient='index')
        meld.index = meld.index.astype(float)
        meld.sort_index(inplace=True)
        meld.set_index([self.data.index.get_level_values(0),self.data.index.get_level_values(1)], inplace=True)
        meld.columns = self.data.columns
        meld.to_csv(outname,sep='\t')


class sva:

    def __init__(self, filename):
        self.data = pd.read_csv(filename,sep='\t').set_index(['Peptide','Protein'])
        self.notdone = True

    def get_tpoints(self):
        tpoints = [i.replace('CT','') for i in self.data.columns.values]
        tpoints = [int(i.split('_')[0]) for i in tpoints]
        self.tpoints = np.asarray(tpoints)

    def prim_cor(self):

        def autocorr(l,shift):
            return dot(l, np.roll(l, shift)) / dot(l, l)

        per = 12
        cors = []
        for row in tqdm(self.data.values):
            ave = []
            for k in set(self.tpoints):
                ave.append((np.mean([row[i] for i, j in enumerate(self.tpoints) if j == k])*1000000))
            cors.append((autocorr(ave,per) - autocorr(ave,(per//2))))
        self.cors = np.asarray(cors)

    def reduce(self,percsub):
        percsub = float(percsub)
        uncor = [(i<(np.percentile(self.cors,percsub))) for i in self.cors]
        self.data_reduced = self.data[uncor]

    def get_res(self):
        res = []
        for row in tqdm(self.data_reduced.values):
            ys = lowess(row, self.tpoints,delta=4)[:,1]
            res.append(row - ys)
        self.res = res

    def get_tks(self):
        pca = RandomizedPCA()
        pca.fit(self.res)
        self.tks = pca.explained_variance_ratio_

    def perm_test(self,nperm):
        nperm = int(nperm)
        def get_res(arr,l):
            res = []
            for row in arr:
                ys = lowess(row, l,delta=4)[:,1]
                res.append(row - ys)
            return res

        def get_tks(arr):
            pca = RandomizedPCA()
            pca.fit(arr)
            return pca.explained_variance_ratio_

        rstar = np.copy(self.res)
        out = np.zeros(len(self.tks))
        for j in tqdm(range(nperm)):
            for i in range(rstar.shape[0]):
                np.random.shuffle(rstar[i,:])
            resstar = get_res(rstar,self.tpoints)
            tkstar = get_tks(resstar)
            #tkstar = get_tks(rstar)
            for m in range(len(self.tks)):
                if tkstar[m] > self.tks[m]:
                    out[m] += 1
        self.sigs = out/nperm

    def eig_reg(self,alpha):
        alpha = float(alpha)
        U, s, V = np.linalg.svd(self.res)
        sig = V.T[:,:len([i for i in itertools.takewhile(lambda x: x < alpha, self.sigs)])]
        pvals = []
        for trend in tqdm(sig.T):
            temp = []
            for row in self.data_reduced.values:
                slope, intercept, r_value, p_value, std_err = linregress(row,trend)
                temp.append(p_value)
            pvals.append(temp)
        self.ps =  pvals

    def subset_svd(self,lam):

        def est_pi_naught(probs_naught,lam):
            return len([i for i in probs_naught if i > lam])/(len(probs_naught)*(1-lam))

        def est_pi_sig(probs_sig,l):
            pi_0 = est_pi_naught(probs_sig,l)
            if pi_0 > 1:
                return 'nan'
            sp = np.sort(probs_sig)
            return sp[int(floor((1-pi_0)*len(probs_sig)))]

        _, _, bt = np.linalg.svd(self.res)
        trends = []
        for j, entry in enumerate(tqdm(self.ps)):
            sub = []
            thresh = est_pi_sig(entry,lam)
            if thresh == 'nan':
                self.ts = trends
            for i in range(len(entry)):
                if entry[i] < thresh:
                    sub.append(self.data_reduced.values[i])
            U, s, V = np.linalg.svd(sub)
            temp = []
            for trend in V:
                _, _, _, p_value, _ = linregress(bt[j],trend)
                temp.append(p_value)
            trends.append(V.T[:,np.argmin(temp)])
        self.ts = trends

    def normalize(self,outname):
        fin_res = np.dot(np.linalg.lstsq(np.asarray(self.ts).T,self.data.values.T)[0].T,np.asarray(self.ts))
        self.svd_norm = self.data.values - fin_res
        self.svd_norm = pd.DataFrame(self.svd_norm,index=self.data.index,columns=self.data.columns)
        self.svd_norm = pd.DataFrame(scale(self.svd_norm.values,axis=1),columns=self.svd_norm.columns,index=self.svd_norm.index)
        self.svd_norm = self.svd_norm.groupby(level='Protein').mean()
        self.svd_norm.index.names = ['#']
        self.svd_norm.to_csv(outname,sep='\t')
