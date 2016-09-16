from limbr import sva
import sys
import getopt

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
            perm = arg
        elif opt in ("-a", "--alpha"):
            a = arg
    print('reading data')
    to_sva = sva(inputfile)
    to_sva.get_tpoints()
    print('calculating primary trend correlations')
    to_sva.prim_cor()
    print('reducing')
    to_sva.reduce(psub)
    print('calculating residuals')
    to_sva.get_res()
    print('calculating explained variance ratios')
    to_sva.get_tks()
    print('permutation testing')
    to_sva.perm_test(perm)
    print('regressing eigentrends')
    to_sva.eig_reg(a)
    print('performing subset SVD')
    l = 0.5
    to_sva.subset_svd(l)
    print('normalizing')
    to_sva.normalize(outputfile)

if __name__ == '__main__':
    main(sys.argv[1:])
