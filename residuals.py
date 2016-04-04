import pandas as pd
import numpy as np

data = pd.read_csv('imputed_peptide_final.txt',sep='\t')


#need tp take list of classes iterate through the list and for the indexes sharing a class average each row and subtract that row from those columns

#then need to do svd on the resulting residuals matrix do significance permutation tests etc... 
