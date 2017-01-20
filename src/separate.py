import numpy as np
import pandas as pd

data = pd.read_csv('./output/actual/imputed_peptide.txt',sep='\t')
data = data.set_index(['Peptide','Protein'])
WT = [i for i in data.columns if 'WT' in i]
C = [i for i in data.columns if 'C' in i]

wt = data[WT]
csp = data[C]

wt.to_csv('./output/actual/wt_for_sva.txt',sep='\t')
csp.to_csv('./output/actual/csp_for_sva.txt',sep='\t')
