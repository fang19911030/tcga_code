# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:24:23 2016

@author: pengcheng
"""

import pandas as pd
import numpy as np
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import MinMaxScaler

"""
get clinical information and count number of M0 and M1
"""
def get_Metastasis(geneID_list,df_clinical):
    num_M0=0
    num_M1=0
    list_Metastasis=[]
    list_Metastasis.append('Metastasis')
    for ID in geneID_list:
        row_number=df_clinical.loc[df_clinical.sampleID==ID].index[0]
        list_Metastasis.append(df_clinical.loc[row_number,'pathologic_M'])
    for i in range(1,len(list_Metastasis)):
        if list_Metastasis[i]!="M0" and isinstance(list_Metastasis[i],str):
            list_Metastasis[i]="M1"
            num_M1+=1
        elif list_Metastasis[i]=="M0":
            num_M0+=1
#        else:
#            list_Metastasis[i]=None
            
    return list_Metastasis, num_M1, num_M0  
            

"""
read Data
"""
path_for_clinical_data="clinical_data"                            
path_for_genomicMatrix="genomicMatrix"
df_clinical = pd.read_table(path_for_clinical_data)
df_genomicMatrix = pd.read_table(path_for_genomicMatrix)

"""
Add metastasis information
"""
geneID_list=list(df_genomicMatrix.columns.values)
geneID_list=geneID_list[1:]
list_M=[]

list_M,num_M1,num_M0=get_Metastasis(geneID_list,df_clinical)
df_genomicMatrix.loc[20531]=list_M

"""
divide train group and test group  and feature normalize
"""
mms=MinMaxScaler()
df_genomicMatrix=df_genomicMatrix.dropna(axis=1)       ##qet rid of sample without metastasis 
X=df_genomicMatrix.iloc[1:20529,1:].values.T           ## one row represent one sample 
y=df_genomicMatrix.iloc[-1,1:].values.T
X_train, X_test, y_train, y_test=train_test_split(X,y,test_size=0.3,random_state=0)
X_train_norm=mms.fit_transform(X_train)
X_test_norm=mms.transform(X_test)


