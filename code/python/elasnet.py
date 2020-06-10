from pathlib import Path
data_folder = Path('C:/Users/langzx/Desktop/github/gene/data')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from functools import reduce
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score,confusion_matrix
from sklearn.metrics import recall_score
from sklearn.model_selection import KFold,cross_val_score
Calu6= pd.read_csv (data_folder/'Calu6.csv')
Calu6['Cell'] = 'Calu6'
#Calu6_M['Signal'] = 1

HCT116= pd.read_csv (data_folder/'HCT116.csv')
HCT116['Cell'] = 'HCT116'
#HCT116_M['Signal'] = 1

LOVO= pd.read_csv (data_folder/'LOVO.csv')
LOVO['Cell'] = 'LOVO'
#LOVO_M['Signal'] = 1

MB231= pd.read_csv (data_folder/'MB231.csv')
MB231['Cell'] = 'MB231'
#MB231_M['Signal'] = 1


MB453= pd.read_csv (data_folder/'MB453.csv')
MB453['Cell'] = 'MB453'
#MB453_M['Signal'] = 0
MB453.head(4)

DLD1= pd.read_csv (data_folder/'DLD1.csv')
DLD1['Cell'] = 'DLD1'
#DLD1['Signal'] = 0
DLD1.head(4)


df_com = [Calu6,HCT116,LOVO,MB231,MB453,DLD1]

dfs = pd.concat(df_com)
dfs.head(4)
dfs.shape


dfs.to_csv(data_folder/"data_concat.csv", index = False)

dd= dfs.reset_index().pivot_table(index = 'Cell', columns = 'ID_REF', values = 'VALUE').fillna(0)
dd.columns
dd.shape

cc= dfs.groupby(['ID_REF']).size().to_frame().reset_index()
cc.columns = ['ID_REF','Size']
repp = cc.loc[cc['Size']>1,]
repp.to_csv(data_folder/'replica23.csv', index = False)


mergedCC = cc.merge(dfs, left_on='ID_REF', right_on='ID_REF', how='inner')
mergedCC.to_csv(data_folder/"replica.csv", index = False)






X = dd
y = np.array([1,1,1,1,0,0])

X_train, X_test, y_train, y_test = \
	train_test_split(X, y, test_size = .2, random_state = 0)


from sklearn.linear_model import LogisticRegression

C = 0.1
l1_ratio = 0.5

clf_en_LR = LogisticRegression(C=C, penalty='elasticnet',
                                   l1_ratio=l1_ratio, tol=0.01,random_state=0,solver= 'saga')

y_pred = clf_en_LR.fit(X_train, y_train).predict(X_test) 


clf = clf_en_LR.fit(X_train, y_train)

clf.coef_


prec = precision_score(y_test, y_pred, average = 'weighted')
recall =recall_score(y_test, y_pred, average = 'weighted')
print("Precision: %s, Recall: %s" %(prec, recall))


coef_table = pd.DataFrame(list(dd.columns)).copy()
coef_table.insert(len(coef_table.columns),"Coefs",clf.coef_.transpose())

coef_table.head(4)

coef_table.columns = ['ID_REF', 'Coefs']
merged = coef_table.merge(mergedCC, left_on='ID_REF', right_on='ID_REF', how='outer').fillna(0)

merged.to_csv(data_folder/"coef_demo_expression.csv", index = False)

# Set up k-fold
k_fold = KFold(n_splits = 4, random_state = 0)

# Evaluate precision and recall for each fold
precision = cross_val_score(
  clf_en_LR, X_train, y_train, cv =k_fold, scoring = 'precision_weighted')
recall = cross_val_score(
  clf_en_LR, X_train, y_train, cv = k_fold, scoring = 'recall_weighted')
print("Precision scores: %s" %(precision)) 
print("Recall scores: %s" %(recall))
# Evaluate precision and recall for each fold

precision = cross_val_score(
  clf, X_train, y_train, cv =k_fold, scoring = 'precision_weighted')
recall = cross_val_score(
  clf, X_train, y_train, cv = k_fold, scoring = 'recall_weighted')
print("Precision scores: %s" %(precision)) 
print("Recall scores: %s" %(recall))

# Iterate over different levels of max depth
for cc in [0.01,0.05,0.1,0.5]:
  # Create and fit model
  clf = LogisticRegression(C=cc, penalty='elasticnet',
                                   l1_ratio=l1_ratio, tol=0.01,random_state=0,solver= 'saga')

  print("Evaluating elastic net = %s" %(cc))
  y_pred = clf.fit(X_train, y_train).predict(X_test) 
  
  # Evaluate confusion matrix, precision, recall
  print("Confusion matrix: ")
  print(confusion_matrix(y_test, y_pred))
  prec = precision_score(y_test, y_pred, average = 'weighted')
  recall = recall_score(y_test,y_pred, average = 'weighted')
  print("Precision: %s, Recall: %s" %(prec, recall))

















dfs['ID_REF'].nunique()
cc = dfs.groupby(['ID_REF']).size().to_frame().reset_index()
cc.columns

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

X = onehot_encoded
X



# Define training and testing
X_train, X_test, y_train, y_test = \
	train_test_split(X, y, test_size = .2, random_state = 0)

from sklearn.linear_model import LogisticRegression

C = 0.05
l1_ratio = 0.5

clf_en_LR = LogisticRegression(C=C, penalty='elasticnet', solver='saga',
                                   l1_ratio=l1_ratio, tol=0.01)

y_pred = clf_en_LR.fit(X_train, y_train).predict(X_test) 

clf = LogisticRegression(random_state=0).fit(X_train, y_train)

clf.coef_
coef_table = pd.DataFrame(list(cc['ID_REF'])).copy()
coef_table.insert(len(coef_table.columns),"Coefs",clf.coef_.transpose())

coef_table.head(4)

