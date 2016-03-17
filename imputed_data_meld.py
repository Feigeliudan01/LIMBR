import pandas as pd
import json

out = {}
for i in range(0,7):
    tempdict = json.loads(open('results_'+str(i)+'.txt').read())
    out.update(tempdict)

out['0']

meld = pd.DataFrame.from_dict(out,orient='index')
meld.index = meld.index.astype(float)
meld.sort_index(inplace=True)
meld.to_csv('imputed_melded.txt',sep='\t')

meld.head(n=10)
