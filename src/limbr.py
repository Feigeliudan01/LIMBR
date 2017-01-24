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

    def __init__(self, filename,design,data_type,blocks=None):
        np.random.seed(4574)
        self.data_type = str(data_type)
        if self.data_type == 'p':
            self.raw_data = pd.read_csv(filename,sep='\t').set_index(['Peptide','Protein'])
        if self.data_type == 'r':
            self.raw_data = pd.read_csv(filename,sep='\t').set_index('#')
        self.designtype = str(design)
        if self.designtype == 'b':
            self.block_design = pickle.load( open( blocks, "rb" ) )
        self.notdone = True

    def pool_normalize(self):
        def gen_norm_dict(l):
            newd = {}
            for i in range(len(l)):
                newd[l[i]] = int(np.ceil((i+1)/5))
            return newd

        def pool_normalize(df,dmap):
            newdf = pd.DataFrame(index=df.index)
            for column in df.columns.values:
                newdf[column] = df[column].div(df['pool_'+'%02d' % dmap[column]],axis='index')
            nonpool = [i for i in newdf.columns if 'pool' not in i]
            newdf = newdf[nonpool]
            return newdf

        def qnorm(df):
            ref = pd.concat([df[col].sort_values().reset_index(drop=True) for col in df], axis=1, ignore_index=True).mean(axis=1).values
            for i in range(0,len(df.columns)):
                df = df.sort_values(df.columns[i])
                df[df.columns[i]] = ref
            return df.sort_index()
        if self.data_type == 'r':
            self.data = qnorm(self.raw_data)
            self.scaler = preprocessing.StandardScaler().fit(self.data.values.T)
            self.data = pd.DataFrame(self.scaler.transform(self.data.values.T).T,columns=self.data.columns,index=self.data.index)
        else:
            norm_map = gen_norm_dict(self.raw_data.columns.values)
            self.data_pnorm = pool_normalize(self.raw_data,norm_map)
            self.data_pnorm = self.data_pnorm.replace([np.inf, -np.inf], np.nan)
            self.data_pnorm = self.data_pnorm.dropna()
            self.data_pnorm = self.data_pnorm.sort_index(axis=1)
            self.data_pnorm = qnorm(self.data_pnorm)
            self.scaler = preprocessing.StandardScaler().fit(self.data_pnorm.values.T)
            self.data = pd.DataFrame(self.scaler.transform(self.data_pnorm.values.T).T,columns=self.data_pnorm.columns,index=self.data_pnorm.index)


    def get_tpoints(self):
        tpoints = [i.replace('CT','') for i in self.data.columns.values]
        tpoints = [int(i.split('_')[0]) for i in tpoints]
        #tpoints = [int(i.split('.')[0]) for i in tpoints]
        self.tpoints = np.asarray(tpoints)

    def prim_cor(self):
        def circ_cor():
            def autocorr(l,shift):
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
            cors = []
            for row in tqdm(self.data.values):
                ys = lowess(row, self.tpoints, it=1)[:,1]
                cors.append(-sum((row - ys)**2))
            self.cors = np.asarray(cors)

        def block_cor():
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
        percsub = float(percsub)
        uncor = [(i<(np.percentile(self.cors,percsub))) for i in self.cors]
        self.data_reduced = self.data[uncor]


    def get_res(self,in_arr):
        def get_l_res(arr):
            res = []
            for row in arr:
                ys = lowess(row, self.tpoints, it=1)[:,1]
                res.append(row - ys)
            return np.array(res)

        def get_b_res(arr):
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

    def set_res(self):
        self.res = self.get_res(self.data_reduced.values)

    def get_tks(self,arr):
        pca = PCA(svd_solver='randomized',random_state=4574)
        pca.fit(arr)
        return pca.explained_variance_ratio_

    def set_tks(self):
        self.tks = self.get_tks(self.res)

    def perm_test(self,nperm,npr):
        mgr = Manager()
        output = mgr.list()
        def single_it(rseed):
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
        alpha = float(alpha)
        U, s, V = np.linalg.svd(self.res)
        #this takewhile might not be working, need to check
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

        def est_pi_naught(probs_naught,lam):
            return len([i for i in probs_naught if i > lam])/(len(probs_naught)*(1-lam))

        def est_pi_sig(probs_sig,l):
            pi_0 = est_pi_naught(probs_sig,l)
            if pi_0 > 1:
                return 'nan'
            sp = np.sort(probs_sig)
            return sp[int(np.floor((1-pi_0)*len(probs_sig)))]

        _, _, bt = np.linalg.svd(self.res)
        trends = []
        for j, entry in enumerate(tqdm(self.ps)):
            sub = []
            thresh = est_pi_sig(entry,lam)
            if thresh == 'nan':
                self.ts = trends
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
        self.ts = trends

    def normalize(self,outname):
        #self.raw_data.mean(axis=1).to_csv(outname.split('.txt')[0]+'_mean.txt',sep='\t')
        #(self.data.std(axis=1)/self.data.mean(axis=1)).to_csv(outname.split('.txt')[0]+'_cv.txt',sep='\t')
        pd.DataFrame(self.ts,columns=self.data.columns).to_csv(outname.split('.txt')[0]+'_trends.txt',sep='\t')
        pd.DataFrame(self.sigs).to_csv(outname.split('.txt')[0]+'_perms.txt',sep='\t')
        pd.DataFrame(self.tks).to_csv(outname.split('.txt')[0]+'_tks.txt',sep='\t')
        fin_res = np.dot(np.linalg.lstsq(np.asarray(self.ts).T,self.data.values.T)[0].T,np.asarray(self.ts))
        #self.svd_norm = self.data.values - fin_res
        self.svd_norm = self.scaler.inverse_transform((self.data.values - fin_res).T).T
        self.svd_norm = pd.DataFrame(self.svd_norm,index=self.data.index,columns=self.data.columns)
        self.svd_norm = self.svd_norm.mul((self.raw_data[self.raw_data.index.isin(self.data.index)].mean(axis=1)/self.svd_norm.mean(axis=1)),axis=0)
        #self.svd_norm = pd.DataFrame(scale(self.svd_norm.values,axis=1),columns=self.svd_norm.columns,index=self.svd_norm.index)
        if self.data_type == 'p':
            self.svd_norm = self.svd_norm.groupby(level='Protein').mean()
        self.svd_norm.index.names = ['#']
        self.svd_norm.to_csv(outname,sep='\t')
