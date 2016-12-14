import numpy as np
import pandas as pd
from sklearn.preprocessing import scale, MinMaxScaler
import pickle

def gen_norm_dict(l):
    newd = {}
    for i in range(len(l)):
        newd[l[i]] = int(np.ceil((i+1)/10)) + 1
    return newd

def pool_normalize(df,dmap):
    newdf = pd.DataFrame(index=df.index)
    for column in df.columns.values:
        if column[0] == 'W' and column[4] != 'o':
            newdf[column] = df[column].div(df['WT_pool_'+'%02d' % dmap[column]],axis='index')
        if column[0] == 'C' and column[4] != 'o':
            newdf[column] = df[column].div(df['C_pool_'+'%02d' % dmap[column]],axis='index')
    return newdf

def qnorm(df):
    ref = pd.concat([df[col].sort_values().reset_index(drop=True) for col in df], axis=1, ignore_index=True).mean(axis=1).values
    for i in range(0,len(df.columns)):
        df = df.sort_values(df.columns[i])
        df[df.columns[i]] = ref
    return df.sort_index()

data = pd.read_csv('./output/actual/imputed_peptide.txt',sep='\t')
data = data.set_index(['Peptide','Protein'])
norm_map = gen_norm_dict(data.columns.values)
data = pool_normalize(data,norm_map)
data = data.replace([np.inf, -np.inf], np.nan)
data = data.dropna()
data = data.sort_index(axis=1)
#data = qnorm(data)

csp_norm = data[data.columns.values[:75]]
wt_norm = data[data.columns.values[75:]]

wt_norm.columns = [i.replace('WT-','CT') for i in wt_norm.columns.values]

wt_norm.columns = [i.replace('-','_') for i in wt_norm.columns.values]

wt_norm = wt_norm[wt_norm.columns.values[3:]]

csp_norm.columns = [i.replace('C-','CT') for i in csp_norm.columns.values]

csp_norm.columns = [i.replace('-','_') for i in csp_norm.columns.values]

csp_norm = csp_norm[csp_norm.columns.values[3:]]

wt_norm = pd.DataFrame(scale(wt_norm.values,axis=1),columns=wt_norm.columns,index=wt_norm.index)

csp_norm = pd.DataFrame(scale(csp_norm.values,axis=1),columns=csp_norm.columns,index=csp_norm.index)

wt_norm.to_csv('./output/actual/wt_for_sva.txt',sep='\t')
csp_norm.to_csv('./output/actual/csp_for_sva.txt',sep='\t')

block_design = [j for i in range(1,12) for j in [i]*3]*2 + [j for i in range(1,3) for j in [i]*3]
pickle.dump(block_design, open( "../output/actual/block_design.p", "wb" ) )

rna =  pd.read_csv('../data/Jen_rnaseq_formatted_raw_counts.txt',sep='\t')
rna = rna.set_index('Transcript')
rna.index.names = ['#']
#rna = qnorm(rna)
rna.to_csv('../output/actual/rna_for_sva.txt',sep='\t')
