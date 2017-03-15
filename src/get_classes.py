import pandas as pd
import sys
import getopt
import optparse

def main(argv):
    parser = optparse.OptionParser()
    parser.add_option('-l','--list', dest='inputfiles', nargs=3 )
    (options, args) = parser.parse_args()
    ch = pd.read_csv(options.inputfiles[1],sep='\t',index_col=0)
    ch = ch[ch.max(axis=1)!=0]
    ch = ch[ch.min(axis=1)!=0]
    okrows = ch.index
    cl = pd.read_csv(options.inputfiles[0],sep='\t',index_col=0)
    cl.ix[okrows].to_csv(options.inputfiles[2],header=False,index=False,sep='\t')


if __name__ == '__main__':
    main(sys.argv[1:])
