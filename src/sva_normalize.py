from limbr import sva
import sys
import getopt

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"h:i:o:s:p:a:d:b:",["help","ifile=","ofile=","sub=","perm=","alpha=","design=","blocks="])
    except getopt.GetoptError:
        print 'residuals.py -i <inputfile> -o <outputfile> -s <subset%> -p <#permutations> -a <alphalevel> -d <designtype> -b <bdesignpath>'
        sys.exit(2)
    b = None
    for opt, arg in opts:
        if opt in ('-h',"--help"):
            print 'residuals.py -i <inputfile> -o <outputfile> -s <subset%> -p <#permutations> -a <alphalevel> -d <designtype> -b <bdesignpath>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-s", "--sub"):
            psub = arg
        elif opt in ("-p", "--perm"):
            perm = arg
        elif opt in ("-a", "--alpha"):
            a = arg
        elif opt in ("-d", "--design"):
            d = arg
        elif opt in ("-b", "--bpath"):
            b = arg
    print('reading data')
    to_sva = sva(inputfile,d,b)
    to_sva.get_tpoints()
    print('calculating primary trend correlations')
    to_sva.prim_cor()
    print('reducing')
    to_sva.reduce(psub)
    print('calculating residuals')
    to_sva.set_res()
    print('calculating explained variance ratios')
    to_sva.set_tks()
    print('permutation testing')
    to_sva.perm_test(perm)
    print('\nregressing eigentrends')
    to_sva.eig_reg(a)
    print('performing subset SVD')
    l = 0.5
    to_sva.subset_svd(l)
    print('normalizing')
    to_sva.normalize(outputfile)

if __name__ == '__main__':
    main(sys.argv[1:])