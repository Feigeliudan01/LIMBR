import pandas as pd

tmtids = pd.read_csv('tmtids.csv')
tmtvac = pd.read_csv('tmt_speedvac_group.csv')
hphdate = pd.read_csv('hphdate.csv')

tmt_merged = pd.merge(tmtids, isodate, left_on='TMT-10 Set', right_on='TMT Set')
tmt_merged = pd.merge(tmt_merged, tmtvac, on='TMT Set')
tmt_merged = pd.merge(tmt_merged, hphdate, on='TMT Set')

tmt_merged = tmt_merged.drop('TMT-10 Set',axis=1)

prepids = pd.read_csv('prepids.csv')
neg20group = pd.read_csv('neg20group.csv')
gilsongroup = pd.read_csv('gilsongroup.csv')
lysegroup = pd.read_csv('lysegroup.csv')
prepdat = pd.read_csv('prep_data.csv')

prep_merged = pd.merge(prepids, neg20group, on='Prep Set')
prep_merged = pd.merge(prep_merged, gilsongroup, on='Prep Set')
prep_merged = pd.merge(prep_merged, lysegroup, on='Prep Set')
prep_merged = pd.merge(prep_merged, prepdat, on='Prep Set')
all_merged = pd.merge(prep_merged, tmt_merged, on='Sample ID')

keep = []
for i in all_merged.index:
    if all_merged.iloc[i]['Sample ID'].split('-')[0] == 'WT':
        keep.append(i)

wt_class = all_merged.iloc[keep]
wt_class['Sample ID'] = wt_class['Sample ID'].map(lambda x: 'CT'+x.split('-')[1]+'_'+x.split('-')[2])

wt_class = wt_class.set_index('Sample ID')
wt_class.to_csv('wt_classes.csv')

keep = []
for i in all_merged.index:
    if all_merged.iloc[i]['Sample ID'].split('-')[0] == 'C':
        keep.append(i)

csp_class = all_merged.iloc[keep]
csp_class['Sample ID'] = csp_class['Sample ID'].map(lambda x: 'CT'+x.split('-')[1]+'_'+x.split('-')[2])

csp_class = csp_class.set_index('Sample ID')
csp_class.to_csv('csp_classes.csv')
