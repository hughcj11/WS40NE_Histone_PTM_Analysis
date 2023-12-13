#Author: Chelsea Hughes
import csv
import re
import numpy as np
from pandas import read_csv
import pandas as pd
from scipy.stats import f_oneway
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from numpy import inf
base_path = "/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/"

#####Generating values for analysis
##Unique hPTMs
#The below code generates the relative abundance value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
sample=[f"Sample {i+1}" for i in range(0,48)]
datastats[sample]= (data.iloc[:,9:57].to_numpy() /(data.iloc[:,58:106])*100).to_numpy()
datastats=datastats.replace(np.nan, 0)
datastats.to_csv(base_path+'/Full_Relative_Abundance.csv', index=False)
#The below code generates the beta value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
sample=[f"Sample {i+1}" for i in range(0,48)]
datastats[sample]= data.iloc[:,9:57].to_numpy() /(data.iloc[:,58:106]+100).to_numpy()
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/Full_beta_values.csv', index=False)


#The below code generates the m value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
sample=[f"Sample {i+1}" for i in range(0,48)]
datastats[sample]=np.log2((data.iloc[:,9:57].to_numpy()/((data.iloc[:,58:106]+100).to_numpy()))/(1-(data.iloc[:,9:57].to_numpy()/((data.iloc[:,58:106]+100).to_numpy()))))
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/Full_M_Values.csv', index=False)

##Global calculations
#The below code generates the global relative abundance value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
sample=[f"Sample {i+1}" for i in range(0,48)]
datastats[sample]= (data.iloc[:,1:49].to_numpy()/(data.iloc[:,50:98])*100).to_numpy()
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/WS40NE Anoxia Total PTM.csv', index=False)







####Statistical analysis
##Unique hPTMs
#The below code performs statistics on the relative abundance data for each unique hPTM (residue,PTM,and histone)

#################Normoxic vs Recovery
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,70:82].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,20:32],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv(base_path+'/DatastatsNormoxicVRecovery.csv', index=False)

#################Normoxic vs 4 days anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["4DAnoxic_Average"]=(data.iloc[:,45:57].sum(axis=1)/data.iloc[:,94:106].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['4DAnoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,44:57],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv(base_path+'/DatastatsNormoxicV4dAnoxia.csv', index=False)

#################Normoxic vs 24h anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["24hAnoxic_Average"]=(data.iloc[:,33:45].sum(axis=1)/data.iloc[:,82:94].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['24hAnoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,32:44],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv(base_path+'/DatastatsNormoxicV24hAnoxia.csv', index=False)








##Global statistics
#The below code performs statistics on the relative abundance data global hPTM changes (independent of residue and histone)

#################Normoxic vs Recovery
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,70:82].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,20:32],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobalNormoxicVRecovery.csv', index=False)

#################Normoxic vs 4 days anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["4DAnoxic_Average"]=(data.iloc[:,45:57].sum(axis=1)/data.iloc[:,94:106].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['4DAnoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,44:57],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobalNormoxicV4DAnoxia.csv', index=False)

#################Normoxic vs 24h anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["24hAnoxic_Average"]=(data.iloc[:,33:45].sum(axis=1)/data.iloc[:,82:94].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['24hAnoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,32:44],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobalNormoxicV24HAnoxia.csv', index=False)




#The below document shows the relative coverage of each modifiable residue (a residue shown as capable of having a PTM). For example, how often K covered by PTMs?
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= pd.DataFrame()
Normoxic_Average=[]
LTAnoxic_Average=[]
STAnoxic_Average=[]
Recovery_Average=[]
GroupedAATotal=data.drop_duplicates('Amino Acid + Position')
GroupedAA=data.groupby(['Amino Acid']).sum().reset_index()
datastats['Amino Acid']=GroupedAA['Amino Acid']
GroupedAATotal_dict={}
for aminoacid, values in zip(GroupedAATotal['Amino Acid'], zip(GroupedAATotal.iloc[:,59:106].values.tolist())):
    values=values[0]
    if aminoacid in GroupedAATotal_dict.keys():
        GroupedAATotal_dict[aminoacid].append(values)
    else:
        GroupedAATotal_dict[aminoacid]=[]
        GroupedAATotal_dict[aminoacid].append(values)
       
for index, row in GroupedAA.iterrows():
    NormoxicNumerator=row.iloc[9:21].sum()
    amino_acid_dict_rows = GroupedAATotal_dict[row['Amino Acid']]
    if len(amino_acid_dict_rows) ==1:
         dict_rows_sum_for_denom = amino_acid_dict_rows[0]
    else:
        dict_rows_sum_for_denom = [sum(i) for i in zip(*amino_acid_dict_rows)]
    NormoxicDenom=sum(dict_rows_sum_for_denom[0:12])
    Normoxic_Average.append((NormoxicNumerator/NormoxicDenom)*100)
    LTAnoxicNumerator=row.iloc[45:57].sum()
    LTAnoxicDenom=sum(dict_rows_sum_for_denom[36:48])
    LTAnoxic_Average.append((LTAnoxicNumerator/LTAnoxicDenom)*100)
    STAnoxicNumerator=row.iloc[33:45].sum()
    STAnoxicDenom=sum(dict_rows_sum_for_denom[24:36])
    STAnoxic_Average.append((STAnoxicNumerator/STAnoxicDenom)*100)
    RecoveryNumerator=row.iloc[21:33].sum()
    RecoveryDenom=sum(dict_rows_sum_for_denom[12:24])
    Recovery_Average.append((RecoveryNumerator/RecoveryDenom)*100)
datastats["Normoxic_Average"]=Normoxic_Average 
datastats["4d_Anoxic_Average"]=LTAnoxic_Average
datastats["24h_Anoxic_Average"]=STAnoxic_Average
datastats["Recovery_Average"]=Recovery_Average   
datastats.to_csv(base_path+'/ResidueCoverage.csv', index=False)



#The below document shows the relative coverage by a specific PTM for each modifiable residue (a residue shown as capable of having a PTM). For example, how often K covered by Ub?
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= pd.DataFrame()
Normoxic_Average=[]
LTAnoxic_Average=[]
STAnoxic_Average=[]
Recovery_Average=[]
GroupedAATotal=data.drop_duplicates('Amino Acid + Position')
GroupedAA=data.groupby(['Amino Acid',"PTM Description"]).sum().reset_index()
datastats['Amino Acid']=GroupedAA['Amino Acid']
datastats["PTM Description"]=GroupedAA["PTM Description"]
GroupedAATotal_dict={}
for aminoacid, values in zip(GroupedAATotal['Amino Acid'], zip(GroupedAATotal.iloc[:,59:106].values.tolist())):
    values=values[0]
    if aminoacid in GroupedAATotal_dict.keys():
        GroupedAATotal_dict[aminoacid].append(values)
    else:
        GroupedAATotal_dict[aminoacid]=[]
        GroupedAATotal_dict[aminoacid].append(values)
       
for index, row in GroupedAA.iterrows():
    NormoxicNumerator=row.iloc[9:21].sum()
    amino_acid_dict_rows = GroupedAATotal_dict[row['Amino Acid']]
    if len(amino_acid_dict_rows) ==1:
         dict_rows_sum_for_denom = amino_acid_dict_rows[0]
    else:
        dict_rows_sum_for_denom = [sum(i) for i in zip(*amino_acid_dict_rows)]
    NormoxicDenom=sum(dict_rows_sum_for_denom[0:12])
    Normoxic_Average.append((NormoxicNumerator/NormoxicDenom)*100)
    LTAnoxicNumerator=row.iloc[45:57].sum()
    LTAnoxicDenom=sum(dict_rows_sum_for_denom[36:48])
    LTAnoxic_Average.append((LTAnoxicNumerator/LTAnoxicDenom)*100)
    STAnoxicNumerator=row.iloc[33:45].sum()
    STAnoxicDenom=sum(dict_rows_sum_for_denom[24:36])
    STAnoxic_Average.append((STAnoxicNumerator/STAnoxicDenom)*100)
    RecoveryNumerator=row.iloc[21:33].sum()
    RecoveryDenom=sum(dict_rows_sum_for_denom[12:24])
    Recovery_Average.append((RecoveryNumerator/RecoveryDenom)*100)
datastats["Normoxic_Average"]=Normoxic_Average 
datastats["4d_Anoxic_Average"]=LTAnoxic_Average
datastats["24h_Anoxic_Average"]=STAnoxic_Average
datastats["Recovery_Average"]=Recovery_Average 
datastats.to_csv(base_path+'/ResidueCoverageByPTM.csv', index=False)
