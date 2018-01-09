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
        """
        takes input data and missingness threshold and initializes imputable object
        """
        self.data = pd.read_csv(filename,sep='\t')
        self.miss = float(missingness)
        self.pats = {}
        self.notdone = True

    def deduplicate(self):
        """
        removes duplicate peptides

        groups rows by peptide, if a peptide appears in more than one row it is removed
        """
        if self.data[self.data.columns.values[1]][0][-2] == "T":
            self.data[self.data.columns.values[1]] = self.data[self.data.columns.values[1]].apply(lambda x: x.split('T')[0])
            self.data = self.data.groupby(['Peptide','Protein']).mean()
        todrop = []
        for name, group in tqdm(self.data.groupby(level='Peptide')):
            if len(group) > 1:
                todrop.append(name)
        self.data = self.data.drop(todrop)

    def drop_missing(self):
        """removes rows which are missing more data than the user specified missingness threshold"""
        self.miss = np.rint(len(self.data.columns)*self.miss)
        self.data = self.data[self.data.isnull().sum(axis=1)<=self.miss]

    def impute(self,outname):
        """
        imputes missing data with KNN

	takes deduplicated data with missing values removed and an output file name

        outputs the imputed data to the given filename
        """
        def match_pat(l,i):
            """
            finds all missingness patterns present in the dataset

            takes a row of data and its index

            if that row has a new missingness pattern, that pattern is added to the list
            whether the missingness pattern is new or not, the index of that row is assigned to the appropriate missingness pattern
            """
            l = "".join(np.isnan(l).astype(int).astype(str))
            if l not in self.pats.keys():
                self.pats[l] = [i]
            else:
                self.pats[l].append(i)

        def get_patterns(arr):
            """
            calls match_pat on all rows of data
            """
            for ind, val in enumerate(arr):
                match_pat(val,ind)

        def sub_imputer(inds,pattern,origarr,comparr):
            """
            single imputation process for a missingness pattern.

            drops columns missing in a given missingness pattern. Then finds nearest neighbors.  Iterates over rows matching missingness pattern, getting indexes of nearest neighbors, averaging nearest neighbrs and replacing 
missing values with corresponding averages.

            takes indexes of rows in missingness pattern, missingness pattern, the original array of data and the array of rows with no missing values

            returns the matrix of rows matching the missingness pattern with the missing values imputed
            """
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
            """
            calls sub_imputer on each missingness pattern and outputs the results to a dict 
            """
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
