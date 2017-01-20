import numpy as np
import pandas as pd

data = pd.read_csv('./output/actual/imputed_peptide.txt',sep='\t')
data = data.set_index(['Peptide','Protein'])
WT = [i for i in data.columns if 'WT' in i]
C = [i for i in data.columns if 'C' in i]

wt = data[WT]
csp = data[C]

wt.columns = [i.replace('WT-','CT') for i in wt.columns.values]
wt.columns = [i.replace('WT_','') for i in wt.columns.values]
wt.columns = [i.replace('-','_') for i in wt.columns.values]

csp.columns = [i.replace('C-','CT') for i in csp.columns.values]
csp.columns = [i.replace('C_','') for i in csp.columns.values]
csp.columns = [i.replace('-','_') for i in csp.columns.values]

wt.to_csv('./output/actual/wt_for_sva.txt',sep='\t')
csp.to_csv('./output/actual/csp_for_sva.txt',sep='\t')
