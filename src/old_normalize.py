from limbr import old_fashioned
import sys
import getopt

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"h:i:o:s:n:p:a:d:e:pp:",["help","ifile=","ofile=","experiment=","pools="])
    except getopt.GetoptError:
        print('residuals.py -i <inputfile> -o <outputfile> -e <experimenttype> -pp <pdesignpath>')
        sys.exit(2)
    b = None
    for opt, arg in opts:
        if opt in ('-h',"--help"):
            print('residuals.py -i <inputfile> -o <outputfile> -e <experimenttype> -pp <pdesignpath>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-e", "--experiment"):
            e = arg
        elif opt in ("-pp", "--pools"):
            pool = arg
    print('reading data')
    to_old = sva(inputfile,e,pp)
    print('pool normalizing')
    to_old.pool_normalize()
    print('normalizing')
    to_old.normalize(outputfile)

if __name__ == '__main__':
    main(sys.argv[1:])
