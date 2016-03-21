import pandas as pd
import numpy as np
def deduplicate(filename):
    data = pd.read_csv(filename,sep='\t')
    if data[data.columns.values[1]][0][-2] == "T":
        data[data.columns.values[1]] = data[data.columns.values[1]].str.replace('T*', '')
    data = data.groupby(['Peptide','Protein']).mean()
    todrop = []
    for name, group in data.groupby(level='Peptide'):
            if len(group) > 1:
                todrop.append(name)
    dedup = data.drop(todrop)
    return dedup
## make this a function which takes raw file name and max missingness and outputs deduplicated file
def drop_missing(df,missing_percentage):
    miss = np.rint(len(df.columns)*missing_percentage)
    datacomp = df[df.isnull().sum(axis=1)<=miss]
    return datacomp
