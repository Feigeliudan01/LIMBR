import itertools
import pandas as pd
import os

def get_btrends():
    temp = pd.read_csv('mb_denoised_block_1_trends.txt',sep='\t',index_col=0)
    out_df = pd.DataFrame(columns=['#']+list([i.split('_')[0] for i in temp.columns]))
    for it in range(1,21):
        trends = pd.read_csv('mb_trends_' + str(it) + '.txt',sep='\t',index_col=0)
        temp = []
        for j in range(1,(len(trends)+1)):
            for i in list(itertools.combinations([0,1,2], j)):
                temp.append(i)
        combs = []
        for i in temp:
            combs.append(np.sum([trends.values[j] for j in list(i)],axis=0))
        out_df = out_df.append(pd.DataFrame(np.insert(combs, 0, it, axis=1),columns=out_df.columns))
    return out_df

def get_bcors():
    out_df = pd.DataFrame(columns=['iteration','cor'])
    for it in range(1,21):
        trends = pd.read_csv('mb_trends_' + str(it) + '.txt',sep='\t',index_col=0)
        temp = []
        temp.append([pearsonr(np.abs(trends.values[0]),np.abs(trends.values[1]))[0]])
        temp.append([pearsonr(np.abs(trends.values[0]),np.abs(trends.values[2]))[0]])
        temp.append([pearsonr(np.abs(trends.values[1]),np.abs(trends.values[2]))[0]])
        out_df = out_df.append(pd.DataFrame(np.insert(temp, 0, it, axis=1),columns=['iteration','cor']))
    return out_df

def get_pats():
    out_df = pd.DataFrame(columns=['iteration','000','001','010','100','011','110','101','111'])
    for it in range(1,21):
        pats = pd.read_csv('mb_simulated_data_key_' + str(it) + '.txt',sep='\t',index_col=0)
        temp = []
        temp.append(len(pats[(pats['trend1']==0) & (pats['trend2']==0) & (pats['trend3']==0)]))
        temp.append(len(pats[(pats['trend1']==0) & (pats['trend2']==0) & (pats['trend3']==1)]))
        temp.append(len(pats[(pats['trend1']==0) & (pats['trend2']==1) & (pats['trend3']==0)]))
        temp.append(len(pats[(pats['trend1']==1) & (pats['trend2']==0) & (pats['trend3']==0)]))
        temp.append(len(pats[(pats['trend1']==0) & (pats['trend2']==1) & (pats['trend3']==1)]))
        temp.append(len(pats[(pats['trend1']==1) & (pats['trend2']==1) & (pats['trend3']==0)]))
        temp.append(len(pats[(pats['trend1']==1) & (pats['trend2']==0) & (pats['trend3']==1)]))
        temp.append(len(pats[(pats['trend1']==1) & (pats['trend2']==1) & (pats['trend3']==1)]))
        out_df = out_df.append(pd.DataFrame(np.asarray([it]+temp).reshape((1, 9)),columns=['iteration','000','001','010','100','011','110','101','111']))
    return out_df


pats = get_pats()

pd.DataFrame(np.asarray([1, 1240, 1268, 1288, 1271, 1250, 1243, 1225, 1215]).reshape((1, 9)),columns=['iteration','000','001','010','100','011','110','101','111'])

btrends = get_btrends()
btrends.to_csv('btrends_for_eJTK.txt',sep='\t',index=False)
btrends.set_index('#',inplace=True)

minps = pd.read_csv('btrends_for_eJTK__jtkout_GammaP.txt',sep='\t')
minps.set_index('ID',inplace=True)
minps['GammaP'].groupby(minps.index).min()

aucs = pd.read_csv('simdata_mb.csv')
aucs.set_index('Iteration',inplace=True)

aut = pd.merge(minps['GammaP'].groupby(minps.index).min().to_frame(),aucs['Circ_auc'].to_frame(),left_index=True,right_index=True)

cors = np.dot(, )

aut = pd.merge(btrends.sum(axis=1).groupby(btrends.index).sum().to_frame(),aucs['Circ_auc'].to_frame(),left_index=True,right_index=True)
aut = pd.merge(test.groupby(test['iteration']).max(),aucs['Circ_auc'].to_frame(),left_index=True,right_index=True)
#correlate bias trends within experiment

plt.scatter(aut.values[:,0],aut.values[:,1])
plt.show()


df = test

temp = []


combs = []
for i in temp:
    combs.append(np.sum([df.values[j] for j in list(i)],axis=0))


for i in range(7):
    plt.plot(np.asarray(combs[i]).T)
    plt.show()
