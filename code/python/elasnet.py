from pathlib import Path
data_folder = Path('C:/Users/langzx/Desktop/web_gene')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from functools import reduce
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split

Calu6_M= pd.read_csv (data_folder/'Calu6_M.csv')
Calu6_M['Cell'] = 'Calu6'
Calu6_M['Signal'] = 1


HCT116_M= pd.read_csv (data_folder/'HCT116_M.csv')
HCT116_M['Cell'] = 'HCT116'
HCT116_M['Signal'] = 1

LOVO_M= pd.read_csv (data_folder/'LOVO_M.csv')
LOVO_M['Cell'] = 'LOVO'
LOVO_M['Signal'] = 1

MB231_M= pd.read_csv (data_folder/'MB231_M.csv')
MB231_M['Cell'] = 'MB231'
MB231_M['Signal'] = 1


MB453_M= pd.read_csv (data_folder/'MB453_M.csv')
MB453_M['Cell'] = 'MB453'
MB453_M['Signal'] = 0
MB453_M.head(4)

DLD1= pd.read_csv (data_folder/'DLD1.csv')
DLD1['Cell'] = 'DLD1'
DLD1['Signal'] = 0
DLD1.head(4)


df_com = [Calu6_M,HCT116_M,LOVO_M,MB231_M,MB453_M]

dfs = pd.concat(df_com)
dfs.head(4)
dfs.shape
dfs['ID_REF'].nunique()
values = dfs['ID_REF']
y = dfs['Signal']
va = dfs['VALUE']
label_encoder = LabelEncoder()
integer_encoded = label_encoder.fit_transform(values)
print(integer_encoded)
# binary encode
onehot_encoder = OneHotEncoder(sparse=False)
integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
onehot_encoded = onehot_encoder.fit_transform(integer_encoded)
onehot_encoded.shape
np.insert(onehot_encoded, 1, va, axis=1)

va*(onehot_encoded)
X = onehot_encoded

# Define training and testing
X_train, X_test, y_train, y_test = \
	train_test_split(X, y, test_size = .2, random_state = 0)

from sklearn.linear_model import LogisticRegression

C = 0.05
l1_ratio = 0.5

clf_en_LR = LogisticRegression(C=C, penalty='elasticnet', solver='saga',
                                   l1_ratio=l1_ratio, tol=0.01)

df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['ID_REF'],
                                            how='outer'), df_com).fillna('NA')

df_merged

mask453 = dat['Cell'] == 'MB-MDA-453'
mask231 = dat['Cell'] == 'MB-MDA-231'
mask116 = dat['Cell'] == 'HCT116'
masklovo = dat['Cell'] == 'LOVO'
dat453 = dat[mask453]
dat231= dat[mask231]
dat116 = dat[mask116]
datlovo = dat[masklovo]
datlovo.head(5)


df_avg = dat.groupby(['Cell','Hugo Symbol'])['ExAC Allelic Fraction'].mean().reset_index()
df_avg

df_pivot = df_avg.pivot(index='Hugo Symbol', columns='Cell', values='ExAC Allelic Fraction')
df_pivot.loc[:,'Row_Total'] = df_pivot.sum(axis=1)
df_pivot.head(5)
df_clear = df_pivot.loc[df_pivot['Row_Total'] != 0]
df_clear.head(5)
df_clear.to_csv(data_folder/'merged_hugo_symbol.csv')

df_pivot.drop('Row_Total',axis = 1, inplace = True)
ax = sns.heatmap(df_pivot)

plt.show()

dfs = [dat453,dat231,dat116, datlovo]
nan_value = 0
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['Hugo Symbol'],
                                            how='outer'), dfs).fillna(nan_value)

df_merged.head(5)
df_merged.columns
df_merged.to_csv(data_folder/'merged_hugo_symbol.csv', index = False)

len(dat['Cell'].unique())
gennames = gene.unique()
len(gennames)

sigdata = dat.loc[dat['Cell'] != 'MB-MDA-453']
sigdata.shape
sigdata.columns
len(sigdata['Hugo Symbol'].unique())


nosigdata = dat.loc[dat['Cell'] == 'MB-MDA-453']
nosigdata.shape
len(nosigdata['Hugo Symbol'].unique())

sigdata.set_index('Hugo Symbol').stack()
df1 = sigdata.groupby('Hugo Symbol')['Cell'].apply(list).reset_index(name='value_sig')

dd1 = sigdata.groupby(['Hugo Symbol','Cell'])
dd1.head(5)

df1.head(5)


df2 = nosigdata.groupby('Hugo Symbol')['Value'].apply(list).reset_index(name='value_others')

df2.head(5)

dfhugo = df1.merge(df2, left_on='Hugo Symbol', right_on='Hugo Symbol', how='outer').fillna(0)
dfhugo.shape
dfhugo.head(6)
dfhugo.to_csv(data_folder/'hugo_symbol.csv', index = False)
dfhugo.shape
len(dfhugo['Hugo Symbol'].unique())
