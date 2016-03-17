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
from Queue import Queue
import math
from itertools import islice
import json

logger = multiprocessing.log_to_stderr(logging.DEBUG)
logger.setLevel(multiprocessing.SUBDEBUG)
logging.basicConfig(level=logging.INFO)
logger2 = logging.getLogger(__name__)

#logging.basicConfig(filename='mp_debug.log',level=logging.DEBUG)

def find_dups(l):
    newl = []
    dups = []
    for entry in l:
        if entry in newl:
            dups.append(entry)
        else:
            newl.append(entry)
    return dups

def drop_dups(df, dups):
    newdf = df.drop(dups)
    return newdf

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

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
        """ The worker function, invoked in a process. 'nums' is a
            list of numbers to factor. The results are placed in
            a dictionary that's pushed to a queue.
        """
        outdict = {}
        n = 0
        for k in keys:
            outdict[k] = sub_imputer(k,patdict[k], origarr,comparr)
            n += 1
            logging.info("%d / %d" % (n,len(keys)))
        with open('results_%d.txt' % out_q, "w") as f:
                   json.dump(outdict, f)
                   #f.write("\n")
    # Each process will get 'chunksize' nums and a queue to put his out
    # dict into
    out_q = Queue()
    chunksize = int(math.ceil(len(patdict.keys()) / float(nprocs)))
    procs = []

    for i in range(nprocs):
        p = multiprocessing.Process(
                target=worker,
                args=(patdict,patdict.keys()[chunksize * i:chunksize * (i + 1)], origarr, comparr,
                      i))
        procs.append(p)
        logging.info("Starting process %d" % i)
        p.start()

    # # Collect all results into a single result dict. We know how many dicts
    # # with results to expect.
    # resultdict = {}
    # logging.info("getting output")
    # for i in range(nprocs):
    #     resultdict.update(out_q.get())
    #     logging.info("got dict %d" % i)
    # # Wait for all worker processes to finish
    # for p in procs:
    #     logging.info("joining dict %d" % i)
    #     p.join()

    return



if __name__ == '__main__':
    # Define the parameters to test
    data = pd.read_csv('Deduplicated_for_imputation.txt', sep = '\t')
    datacomp = data[data.isnull().sum(axis=1)<=64]
    datavals = datacomp.drop(data.columns[:3],axis=1).values
    comparr = datavals[~np.isnan(datavals).any(axis=1)]
    patterns = get_patterns(datavals)
    revpatterns = {i:x for x,y in patterns.iteritems() for i in y}

    mp_imputer(revpatterns, datavals, comparr, 7)

# data={}
#
# f = open('results.txt')
# for line in iter(f):
#     data.update(json.loads(line))
# f.close()
#
#
# with open('results.txt', "r") as f:
#     while line:
#         line = f.readline()
#         data.update(json.loads(line))
# len(data.keys())
# len(data['51478'])
