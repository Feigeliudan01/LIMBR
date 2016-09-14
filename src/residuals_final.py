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

class sva:

    def __init__(self, filename, percsub, nperm, alpha):
        self.data = pd.read_csv(filename,sep='\t').set_index(['Peptide','Protein'])
        self.sub = float(percsub)
        self.perm = int(nperm)
        self.alpha = float(alpha)
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

    def reduce(self):
        uncor = [(i<(np.percentile(self.cors,self.sub))) for i in self.cors]
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

    def perm_test(self):

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
        for j in tqdm(range(self.perm)):
            for i in range(rstar.shape[0]):
                np.random.shuffle(rstar[i,:])
            resstar = get_res(rstar,self.tpoints)
            tkstar = get_tks(resstar)
            #tkstar = get_tks(rstar)
            for m in range(len(self.tks)):
                if tkstar[m] > self.tks[m]:
                    out[m] += 1
        self.sigs = out/self.perm

    def eig_reg(self):
        U, s, V = np.linalg.svd(self.res)
        sig = V.T[:,:len([i for i in itertools.takewhile(lambda x: x < self.alpha, self.sigs)])]
        pvals = []
        for trend in tqdm(sig.T):
            temp = []
            for row in self.data_reduced.values:
                slope, intercept, r_value, p_value, std_err = linregress(row,trend)
                temp.append(p_value)
            pvals.append(temp)
        self.ps =  pvals

    def subset_svd(self):

        def est_pi_naught(probs_naught,lam):
            return len([i for i in probs_naught if i > lam])/(len(probs_naught)*(1-lam))

        def est_pi_sig(probs_sig,l):
            pi_0 = est_pi_naught(probs_sig,l)

            if pi_0 > 1:
                return 'nan'

            sp = np.sort(probs_sig)
            return sp[int(floor((1-pi_0)*len(probs_sig)))]

        lam = 0.5
        _, _, bt = np.linalg.svd(self.res)
        trends = []
        for j, entry in tqdm(enumerate(self.ps)):
            sub = []
            thresh = est_pi_sig(entry,lam)
            if thresh == 'nan':
                return trends
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

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"h:i:o:s:p:a:",["help","ifile=","ofile=","sub=","perm=","alpha="])
    except getopt.GetoptError:
        print 'residuals.py -i <inputfile> -o <outputfile> -s <subset%> -p <#permutations> -a <alphalevel>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h',"--help"):
            print 'residuals.py -i <inputfile> -o <outputfile> -s <subset%> -p <#permutations> -a <alphalevel>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-s", "--sub"):
            psub = arg
        elif opt in ("-p", "--perm"):
            nperm = arg
        elif opt in ("-a", "--alpha"):
            alpha = arg
    print('reading data')
    to_sva = sva(inputfile,psub,nperm,alpha)
    to_sva.get_tpoints()
    print('calculating primary trend correlations')
    to_sva.prim_cor()
    print('reducing')
    to_sva.reduce()
    print('calculating residuals')
    to_sva.get_res()
    print('calculating explained variance ratios')
    to_sva.get_tks()
    print('permutation testing')
    to_sva.perm_test()
    print('regressing eigentrends')
    to_sva.eig_reg()
    print('performing subset SVD')
    to_sva.subset_svd()
    print('normalizing')
    to_sva.normalize(outputfile)

if __name__ == '__main__':
    main(sys.argv[1:])
