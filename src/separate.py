import numpy as np
import pandas as pd
import pickle

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

block_design = [j for i in range(1,12) for j in [i]*3]*2 + [j for i in range(1,3) for j in [i]*3]
pickle.dump(block_design, open( "./output/actual/block_design.p", "wb" ) )

def qnorm(df):
    ref = pd.concat([df[col].sort_values().reset_index(drop=True) for col in df], axis=1, ignore_index=True).mean(axis=1).values
    for i in range(0,len(df.columns)):
        df = df.sort_values(df.columns[i])
        df[df.columns[i]] = ref
    return df.sort_index()

def gen_norm_dict(l):
    newd = {}
    for i in range(len(l)):
        newd[l[i]] = int(np.ceil((i+1)/5))
    return newd

norm_map =

pickle.dump(gen_norm_dict(wt.columns.values), open( "./output/actual/wt_pool_design.p", "wb" ) )
pickle.dump(gen_norm_dict(csp.columns.values), open( "./output/actual/csp_pool_design.p", "wb" ) )

rna =  pd.read_csv('./data/Jen_rnaseq_formatted_raw_counts.txt',sep='\t')
rna = rna.set_index('Transcript')
rna.index.names = ['#']
rna = rna[rna.sum(axis=1)>0]
rna.to_csv('./output/actual/rna_for_sva.txt',sep='\t')
