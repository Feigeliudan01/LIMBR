###GET PATTERNS
# drop all rows without any nas
# loop through array applying isnan to each row and storing in new array
# loop through new array adding row as key in dict and index of row as value if not already in dict
# otherwise append index to value of pre-existing key

##IMPUTE
# return to original array, for each key in dict, drop all columns for which isnan was true
# then drop all rows still containing nans
# then do sklearn knn and get neighbor matrix
# for each value in this key, get indexes of nearest neighbors (1s in neighbors matrix), average those rows
# iterate through entries in row and replace indexes where isnan is true with value from average
# broadcast row to final copy of array (x[0] = np.arange(10))
#return final copy

import numpy as np
from sklearn.neighbors import NearestNeighbors
import pandas as pd
import multiprocessing, logging
#from Queue import Queue
import math
#from itertools import islice
import json
import os
import sys
import getopt
import time
from ctypes import c_int
import timeit

logger = multiprocessing.log_to_stderr(logging.DEBUG)
logger.setLevel(multiprocessing.SUBDEBUG)
logging.basicConfig(level=logging.INFO)
logger2 = logging.getLogger(__name__)
logging.basicConfig(filename='mp_debug.log',level=logging.DEBUG)

class imputable:

    def __init__(self, filename, missingness, nprocs):
        self.data = pd.read_csv(filename,sep='\t')
        self.miss = float(missingness)
        self.nprocs = int(nprocs)
        self.notdone = True

    def deduplicate(self):
        if self.data[self.data.columns.values[1]][0][-2] == "T":
            self.data[self.data.columns.values[1]] = self.data[self.data.columns.values[1]].apply(lambda x: x.split('T')[0])
            self.data = self.data.groupby(['Peptide','Protein']).mean()
        todrop = []
        for name, group in self.data.groupby(level='Peptide'):
            if len(group) > 1:
                todrop.append(name)
        self.data = self.data.drop(todrop)

    def drop_missing(self):
        self.miss = np.rint(len(self.data.columns)*self.miss)
        self.data = self.data[self.data.isnull().sum(axis=1)<=self.miss]

    def impute(self):
        def match_pat(l,d,i):
            l = "".join(np.isnan(l).astype(int).astype(str))
            if l not in d.keys():
                d[l] = [i]
            else:
                d[l].append(i)

        def get_patterns(arr):
            patts = {}
            for ind, val in enumerate(arr):
                match_pat(val,patts,ind)
            return patts

        def sub_imputer(key,pattern,origarr,comparr):
            newarr = comparr[:,~np.array(list(pattern)).astype(bool)]
            nbrs = NearestNeighbors(n_neighbors=10).fit(newarr)
            indexes = nbrs.kneighbors([origarr[key,~np.array(list(pattern)).astype(bool)]],return_distance=False)
            means = np.mean(comparr[indexes[0][1:]], axis=0)
            outl = []
            for ind, v in enumerate(origarr[key]):
                if not np.isnan(v):
                    outl.append(v)
                else:
                    outl.append(means[ind])
            logging.debug("imputed %d" % key)
            return outl

        def mp_imputer(patdict, origarr, comparr, nprocs):
            def worker(patdict, keys, origarr, comparr, out_q):
                """ The worker function, invoked in a process.
                    it takes a
                """
                outdict = {}
                n = 0
                for k in keys:
                    outdict[k] = sub_imputer(k,patdict[k], origarr,comparr)
                    n += 1
                    logging.info("%d / %d" % (n,len(keys)))
                with open('results_%d.txt' % out_q, "w+") as f:
                           json.dump(outdict, f)
                with counter_lock:
                    counter.value += 1
            chunksize = int(math.ceil(len(patdict.keys()) / float(nprocs)))
            self.procs = []

            for i in range(nprocs):
                p = multiprocessing.Process(target=worker,args=(patdict,patdict.keys()[chunksize * i:chunksize * (i + 1)], origarr, comparr,i))
                self.procs.append(p)
                logging.info("Starting process %d" % i)
                p.start()
            return

        datavals = self.data.values
        comparr = datavals[~np.isnan(datavals).any(axis=1)]
        patterns = get_patterns(datavals)
        revpatterns = {i:x for x,y in patterns.iteritems() for i in y}
        mp_imputer(revpatterns, datavals, comparr, self.nprocs)

    def check(self):
        while self.notdone:
            time.sleep(60)
            if counter.value == self.nprocs:
                self.notdone = False
                for p in self.procs:
                    p.join()


    def meld(self,outname):
        out = {}
        for i in range(0,self.nprocs):
            tempdict = json.loads(open('results_'+str(i)+'.txt').read())
            os.remove('results_'+str(i)+'.txt')
            out.update(tempdict)
        meld = pd.DataFrame.from_dict(out,orient='index')
        meld.index = meld.index.astype(float)
        meld.sort_index(inplace=True)
        meld.set_index([self.data.index.get_level_values(0),self.data.index.get_level_values(1)], inplace=True)
        meld.columns = self.data.columns
        meld.to_csv(outname,sep='\t')


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:m:p:",["ifile=","ofile=","miss=","proc="])
    except getopt.GetoptError:
        print 'knn_impute.py -i <inputfile> -o <outputfile> -m <missing%> -p <#processors>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'knn_impute.py -i <inputfile> -o <outputfile> -m <missing%> -p <#processors>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-m", "--miss"):
            pmiss = arg
        elif opt in ("-p", "--proc"):
            nprocs = arg
    to_impute = imputable(inputfile,pmiss,nprocs)
    to_impute.deduplicate()
    to_impute.drop_missing()
    to_impute.impute()
    to_impute.check()
    to_impute.meld(outputfile)

if __name__ == '__main__':
    start = timeit.timeit()
    counter = multiprocessing.Value(c_int)  # defaults to 0
    counter_lock = multiprocessing.Lock()
    main(sys.argv[1:])
    end = timeit.timeit()
    print end - start
