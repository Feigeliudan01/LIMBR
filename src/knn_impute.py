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
from limbr import imputable
#import multiprocessing, logging
import os
import sys
import getopt

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:m:p:",["ifile=","ofile=","miss=","proc="])
    except getopt.GetoptError:
        print('knn_impute.py -i <inputfile> -o <outputfile> -m <missing%> -p <#processors>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('knn_impute.py -i <inputfile> -o <outputfile> -m <missing%> -p <#processors>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-m", "--miss"):
            pmiss = arg
        elif opt in ("-p", "--proc"):
            nprocs = arg
    print('Reading Data')
    to_impute = imputable(inputfile,pmiss)
    print('Deduplicating')
    to_impute.deduplicate()
    print('Dropping Rows Over Missing Value Threshold')
    to_impute.drop_missing()
    print('Imputing Missing Values')
    to_impute.impute(outputfile)
    #to_impute.check()
    #to_impute.meld(outputfile)

if __name__ == '__main__':
    main(sys.argv[1:])
