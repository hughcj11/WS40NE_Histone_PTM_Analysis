#Author: Chelsea Hughes
import csv
import re
import numpy as np
from pandas import read_csv
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
#data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
#datastats= data.iloc[:,0:8]
#datastats["Normoxic_Average"]=data.iloc[:,8:14].mean(axis=1)
#datastats["Anoxic_Average"]=data.iloc[:,14:20].mean(axis=1)
#datastats['log2FC'] = np.log2(datastats['Anoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
# datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:14],data.iloc[:,14:20],axis=1)
#datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data.iloc[:,8:14], data.iloc[:,14:20],axis=1)
#datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
#datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]

# Write DataFrame to CSV
#datastats.to_csv(base_path+'/Datastats.csv', index=False)

##Global statistics
#The below code performs statistics on the relative abundance data global hPTM changes (independent of residue and histone)
#data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Embryo Anoxia Total PTM.csv",header=0)
#datastats= data.iloc[:,0:1]
#datastats["Normoxic_Average"]=data.iloc[:,1:7].mean(axis=1)
#datastats["Anoxic_Average"]=data.iloc[:,7:13].mean(axis=1)
#datastats['log2FC'] = np.log2(datastats['Anoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
# datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:14],data.iloc[:,14:20],axis=1)
#datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data.iloc[:,1:7], data.iloc[:,7:13],axis=1)
#datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
#datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
#datastats.to_csv(base_path+'/DatastatsGlobal.csv', index=False)