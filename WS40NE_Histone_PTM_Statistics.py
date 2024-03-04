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
sample=[f"Sample {i+1}" for i in range(0,46)]
#sample=[f"Sample {i+1}" for i in range(0,48)]
datastats[sample]= (data.iloc[:,9:55].to_numpy() /(data.iloc[:,56:104])*100).to_numpy()
#datastats[sample]= (data.iloc[:,9:57].to_numpy() /(data.iloc[:,58:106])*100).to_numpy()
datastats=datastats.replace(np.nan, 0)
datastats.to_csv(base_path+'/Full_Relative_Abundance.csv', index=False)
#The below code generates the beta value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
sample=[f"Sample {i+1}" for i in range(0,46)]
#sample=[f"Sample {i+1}" for i in range(0,48)]
datastats[sample]= (data.iloc[:,9:55].to_numpy() /(data.iloc[:,56:104])*100).to_numpy()
#datastats[sample]= data.iloc[:,9:57].to_numpy() /(data.iloc[:,58:106]+100).to_numpy()
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/Full_beta_values.csv', index=False)


#The below code generates the m value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
sample=[f"Sample {i+1}" for i in range(0,46)]
#sample=[f"Sample {i+1}" for i in range(0,48)]
datastats[sample]=np.log2((data.iloc[:,9:55].to_numpy()/((data.iloc[:,56:104]+100).to_numpy()))/(1-(data.iloc[:,9:55].to_numpy()/((data.iloc[:,56:104]+100).to_numpy()))))
#datastats[sample]=np.log2((data.iloc[:,9:57].to_numpy()/((data.iloc[:,58:106]+100).to_numpy()))/(1-(data.iloc[:,9:57].to_numpy()/((data.iloc[:,58:106]+100).to_numpy()))))
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/Full_M_Values.csv', index=False)

##Global calculations
#The below code generates the global relative abundance value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
#sample=[f"Sample {i+1}" for i in range(0,48)]
sample=[f"Sample {i+1}" for i in range(0,46)]
#datastats[sample]= (data.iloc[:,1:49].to_numpy()/(data.iloc[:,50:98])*100).to_numpy()
datastats[sample]= (data.iloc[:,1:47].to_numpy()/(data.iloc[:,48:94])*100).to_numpy()
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/WS40NE Anoxia Total PTM.csv', index=False)







####Statistical analysis
##Unique hPTMs
#The below code performs statistics on the relative abundance data for each unique hPTM (residue,PTM,and histone)



#################Normoxic vs Recovery
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,56:68].sum(axis=1))*100
#datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,68:80].sum(axis=1))*100
#datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,70:82].sum(axis=1))*100
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
#datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,56:68].sum(axis=1))*100
#datastats["4DAnoxic_Average"]=(data.iloc[:,45:57].sum(axis=1)/data.iloc[:,94:106].sum(axis=1))*100
datastats["4DAnoxic_Average"]=(data.iloc[:,44:55].sum(axis=1)/data.iloc[:,91:102].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['4DAnoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)

datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,43:54],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv(base_path+'/DatastatsNormoxicV4dAnoxia.csv', index=False)

#################Normoxic vs 24h anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
#datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,56:68].sum(axis=1))*100
#datastats["24hAnoxic_Average"]=(data.iloc[:,33:45].sum(axis=1)/data.iloc[:,82:94].sum(axis=1))*100
datastats["24hAnoxic_Average"]=(data.iloc[:,33:44].sum(axis=1)/data.iloc[:,80:91].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['24hAnoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,32:43],axis=1)
#datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,32:44],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv(base_path+'/DatastatsNormoxicV24hAnoxia.csv', index=False)


#################4Danoxia vs 24h anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
#datastats["4DAnoxic_Average"]=(data.iloc[:,45:57].sum(axis=1)/data.iloc[:,94:106].sum(axis=1))*100
datastats["4DAnoxic_Average"]=(data.iloc[:,44:55].sum(axis=1)/data.iloc[:,91:102].sum(axis=1))*100
#datastats["24hAnoxic_Average"]=(data.iloc[:,33:45].sum(axis=1)/data.iloc[:,82:94].sum(axis=1))*100
datastats["24hAnoxic_Average"]=(data.iloc[:,33:44].sum(axis=1)/data.iloc[:,80:91].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['4DAnoxic_Average'] / (datastats['24hAnoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
#datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:20], data2.iloc[:,32:44],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,43:54], data2.iloc[:,32:43],axis=1)
data2.iloc[:,32:43]
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv(base_path+'/Datastats4DV24hAnoxia.csv', index=False)

#################4Danoxia vs Recovery
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
#datastats["4DAnoxic_Average"]=(data.iloc[:,45:57].sum(axis=1)/data.iloc[:,94:106].sum(axis=1))*100
datastats["4DAnoxic_Average"]=(data.iloc[:,44:55].sum(axis=1)/data.iloc[:,91:102].sum(axis=1))*100
#datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,70:82].sum(axis=1))*100
datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,68:80].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery_Average'] / (datastats['4DAnoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,20:32], data2.iloc[:,43:54],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv(base_path+'/DatastatsRecoveryV4DAnoxia.csv', index=False)


#################24hanoxia vs Recovery
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
#datastats["24hAnoxic_Average"]=(data.iloc[:,33:45].sum(axis=1)/data.iloc[:,82:94].sum(axis=1))*100
datastats["24hAnoxic_Average"]=(data.iloc[:,33:44].sum(axis=1)/data.iloc[:,80:91].sum(axis=1))*100
#datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,70:82].sum(axis=1))*100
datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,68:80].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery_Average'] / (datastats['24hAnoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,20:32], data2.iloc[:,32:43],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv(base_path+'/DatastatsRecoveryV24hAnoxia.csv', index=False)







##Global statistics
#The below code performs statistics on the relative abundance data global hPTM changes (independent of residue and histone)

#################Normoxic vs Recovery
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
#datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["Normoxic_Average"]=(data.iloc[:,1:13].sum(axis=1)/data.iloc[:,48:60].sum(axis=1))*100
#datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,70:82].sum(axis=1))*100
datastats["Recovery_Average"]=(data.iloc[:,13:25].sum(axis=1)/data.iloc[:,60:72].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,1:13], data2.iloc[:,13:25],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobalNormoxicVRecovery.csv', index=False)

#################Normoxic vs 4 days anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
#datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["Normoxic_Average"]=(data.iloc[:,1:13].sum(axis=1)/data.iloc[:,48:60].sum(axis=1))*100
#datastats["4DAnoxic_Average"]=(data.iloc[:,45:57].sum(axis=1)/data.iloc[:,94:106].sum(axis=1))*100
datastats["4DAnoxic_Average"]=(data.iloc[:,36:47].sum(axis=1)/data.iloc[:,83:94].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['4DAnoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,1:13], data2.iloc[:,36:47],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobalNormoxicV4DAnoxia.csv', index=False)

#################Normoxic vs 24h anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
#datastats["Normoxic_Average"]=(data.iloc[:,9:21].sum(axis=1)/data.iloc[:,58:70].sum(axis=1))*100
datastats["Normoxic_Average"]=(data.iloc[:,1:13].sum(axis=1)/data.iloc[:,48:60].sum(axis=1))*100
#datastats["24hAnoxic_Average"]=(data.iloc[:,33:45].sum(axis=1)/data.iloc[:,82:94].sum(axis=1))*100
datastats["24hAnoxic_Average"]=(data.iloc[:,25:36].sum(axis=1)/data.iloc[:,72:83].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['24hAnoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,1:13], data2.iloc[:,25:36],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobalNormoxicV24HAnoxia.csv', index=False)


#################4Danoxia vs 24h anoxia
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
#datastats["4DAnoxic_Average"]=(data.iloc[:,45:57].sum(axis=1)/data.iloc[:,94:106].sum(axis=1))*100
datastats["4DAnoxic_Average"]=(data.iloc[:,36:47].sum(axis=1)/data.iloc[:,83:94].sum(axis=1))*100
#datastats["24hAnoxic_Average"]=(data.iloc[:,33:45].sum(axis=1)/data.iloc[:,82:94].sum(axis=1))*100
datastats["24hAnoxic_Average"]=(data.iloc[:,25:36].sum(axis=1)/data.iloc[:,72:83].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['4DAnoxic_Average'] / (datastats['24hAnoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,25:36], data2.iloc[:,36:47],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobal4DV24hAnoxia.csv', index=False)

#################4Danoxia vs Recovery
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
#datastats["4DAnoxic_Average"]=(data.iloc[:,45:57].sum(axis=1)/data.iloc[:,94:106].sum(axis=1))*100
datastats["4DAnoxic_Average"]=(data.iloc[:,36:47].sum(axis=1)/data.iloc[:,83:94].sum(axis=1))*100
#datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,70:82].sum(axis=1))*100
datastats["Recovery_Average"]=(data.iloc[:,13:25].sum(axis=1)/data.iloc[:,60:72].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery_Average'] / (datastats['4DAnoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,36:47], data2.iloc[:,13:25],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobalRecoveryV4DAnoxia.csv', index=False)

#################24hanoxia vs Recovery
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/WS40NE Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/WS40NE Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
#datastats["24hAnoxic_Average"]=(data.iloc[:,33:45].sum(axis=1)/data.iloc[:,82:94].sum(axis=1))*100
datastats["24hAnoxic_Average"]=(data.iloc[:,25:36].sum(axis=1)/data.iloc[:,72:83].sum(axis=1))*100
#datastats["Recovery_Average"]=(data.iloc[:,21:33].sum(axis=1)/data.iloc[:,70:82].sum(axis=1))*100
datastats["Recovery_Average"]=(data.iloc[:,13:25].sum(axis=1)/data.iloc[:,60:72].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery_Average'] / (datastats['24hAnoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
#datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:20],data.iloc[:,20:32],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,25:36], data2.iloc[:,13:25],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobalRecoveryV24hAnoxia.csv', index=False)



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
for aminoacid, values in zip(GroupedAATotal['Amino Acid'], zip(GroupedAATotal.iloc[:,56:106].values.tolist())):
    values=values[0]
    if aminoacid in GroupedAATotal_dict.keys():
        GroupedAATotal_dict[aminoacid].append(values)
    else:
        GroupedAATotal_dict[aminoacid]=[]
        GroupedAATotal_dict[aminoacid].append(values)
       
for index, row in GroupedAA.iterrows():
    amino_acid_dict_rows = GroupedAATotal_dict[row['Amino Acid']]
    if len(amino_acid_dict_rows) ==1:
         dict_rows_sum_for_denom = amino_acid_dict_rows[0]
    else:
        dict_rows_sum_for_denom = [sum(i) for i in zip(*amino_acid_dict_rows)]
    NormoxicNumerator=row.iloc[9:21].sum()
    NormoxicDenom=sum(dict_rows_sum_for_denom[0:12])
    Normoxic_Average.append((NormoxicNumerator/NormoxicDenom)*100)
    LTAnoxicNumerator=row.iloc[44:55].sum()
    LTAnoxicDenom=sum(dict_rows_sum_for_denom[35:46])
    LTAnoxic_Average.append((LTAnoxicNumerator/LTAnoxicDenom)*100)
    STAnoxicNumerator=row.iloc[33:44].sum()
    STAnoxicDenom=sum(dict_rows_sum_for_denom[24:35])
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
for aminoacid, values in zip(GroupedAATotal['Amino Acid'], zip(GroupedAATotal.iloc[:,56:106].values.tolist())):
    values=values[0]
    if aminoacid in GroupedAATotal_dict.keys():
        GroupedAATotal_dict[aminoacid].append(values)
    else:
        GroupedAATotal_dict[aminoacid]=[]
        GroupedAATotal_dict[aminoacid].append(values)
       
for index, row in GroupedAA.iterrows():
    amino_acid_dict_rows = GroupedAATotal_dict[row['Amino Acid']]
    if len(amino_acid_dict_rows) ==1:
         dict_rows_sum_for_denom = amino_acid_dict_rows[0]
    else:
        dict_rows_sum_for_denom = [sum(i) for i in zip(*amino_acid_dict_rows)]
    NormoxicNumerator=row.iloc[9:21].sum()
    NormoxicDenom=sum(dict_rows_sum_for_denom[0:12])
    Normoxic_Average.append((NormoxicNumerator/NormoxicDenom)*100)
    LTAnoxicNumerator=row.iloc[44:55].sum()
    LTAnoxicDenom=sum(dict_rows_sum_for_denom[35:46])
    LTAnoxic_Average.append((LTAnoxicNumerator/LTAnoxicDenom)*100)
    STAnoxicNumerator=row.iloc[33:44].sum()
    STAnoxicDenom=sum(dict_rows_sum_for_denom[24:35])
    STAnoxic_Average.append((STAnoxicNumerator/STAnoxicDenom)*100)
    RecoveryNumerator=row.iloc[21:33].sum()
    RecoveryDenom=sum(dict_rows_sum_for_denom[12:24])
    Recovery_Average.append((RecoveryNumerator/RecoveryDenom)*100)
datastats["Normoxic_Average"]=Normoxic_Average 
datastats["4d_Anoxic_Average"]=LTAnoxic_Average
datastats["24h_Anoxic_Average"]=STAnoxic_Average
datastats["Recovery_Average"]=Recovery_Average 
datastats.to_csv(base_path+'/ResidueCoverageByPTM.csv', index=False)

#Code to find all sig hptms an compile them in one document
DataRecv24= read_csv(base_path+'/DatastatsRecoveryV24hAnoxia.csv') #
DataRecv4d= read_csv(base_path+'/DatastatsRecoveryV4DAnoxia.csv') #
DataNormv24h= read_csv(base_path+'/DatastatsNormoxicV24hAnoxia.csv')#
DataNormv4d= read_csv(base_path+'/DatastatsNormoxicV4dAnoxia.csv')#
DataNormvRec= read_csv(base_path+'/DatastatsNormoxicVRecovery.csv')
Data24hv4d= read_csv(base_path+'/Datastats4DV24hAnoxia.csv') #

DataRecv24Filtered= DataRecv24[DataRecv24['corrected_p_values']<.05]
DataRecv24Filtered.rename(columns={"corrected_p_values": "Recovery vs 24hr Anoxia corrected p-value"}, inplace=True) 
DataRecv24Filtered= DataRecv24Filtered.drop(columns=['Recovery_Average','24hAnoxic_Average', 'T-test_t_statistic', 'T-test_p_value', 'log2FC'])

Data24hv4dFiltered= Data24hv4d[Data24hv4d['corrected_p_values']<.05]
Data24hv4dFiltered.rename(columns={"corrected_p_values": "4d Anoxia vs 24hr Anoxia corrected p-value"}, inplace=True) 
Data24hv4dFiltered= Data24hv4dFiltered.drop(columns=['4DAnoxic_Average','24hAnoxic_Average', 'T-test_t_statistic', 'T-test_p_value', 'log2FC'])

DataNormv24hFiltered= DataNormv24h[DataNormv24h['corrected_p_values']<.05]
DataNormv24hFiltered.rename(columns={"corrected_p_values": "Normoxia vs 24hr Anoxia corrected p-value"}, inplace=True) 
DataNormv24hFiltered= DataNormv24hFiltered.drop(columns=['Normoxic_Average','24hAnoxic_Average', 'T-test_t_statistic', 'T-test_p_value', 'log2FC'])

DataNormv4dFiltered= DataNormv4d[DataNormv4d['corrected_p_values']<.05]
DataNormv4dFiltered.rename(columns={"corrected_p_values": "Normoxia vs 4d Anoxia corrected p-value"}, inplace=True) 
DataNormv4dFiltered= DataNormv4dFiltered.drop(columns=['Normoxic_Average','4DAnoxic_Average', 'T-test_t_statistic', 'T-test_p_value', 'log2FC'])

DataNormvRecFiltered= DataNormvRec[DataNormvRec['corrected_p_values']<.05]
DataNormvRecFiltered.rename(columns={"corrected_p_values": "Normoxia vs Recovery corrected p-value"}, inplace=True) 
DataNormvRecFiltered= DataNormvRecFiltered.drop(columns=['Normoxic_Average','Recovery_Average', 'T-test_t_statistic', 'T-test_p_value', 'log2FC'])

DataRecv4dFiltered= DataRecv4d[DataRecv4d['corrected_p_values']<.05]
DataRecv4dFiltered.rename(columns={"corrected_p_values": "4d Anoxia vs Recovery corrected p-value"},inplace=True) 
DataRecv4dFiltered= DataRecv4dFiltered.drop(columns=['4DAnoxic_Average','Recovery_Average', 'T-test_t_statistic', 'T-test_p_value', 'log2FC'])
merged = pd.merge(DataRecv4dFiltered, DataRecv24Filtered, on=['hPTM_ID', 'PTM Description', "Protein Accession", "Protein Description", "Position", "Amino Acid", "Amino Acid + Position", "Unimod"], how='outer')
merged = pd.merge(merged, Data24hv4dFiltered, on=['hPTM_ID', 'PTM Description', "Protein Accession", "Protein Description", "Position", "Amino Acid", "Amino Acid + Position", "Unimod"], how='outer') 
merged = pd.merge(merged, DataNormv24hFiltered, on=['hPTM_ID', 'PTM Description', "Protein Accession", "Protein Description", "Position", "Amino Acid", "Amino Acid + Position", "Unimod"], how='outer') 
merged = pd.merge(merged, DataNormv4dFiltered, on=['hPTM_ID', 'PTM Description', "Protein Accession", "Protein Description", "Position", "Amino Acid", "Amino Acid + Position", "Unimod"], how='outer') 
merged = pd.merge(merged, DataNormvRecFiltered, on=['hPTM_ID', 'PTM Description', "Protein Accession", "Protein Description", "Position", "Amino Acid", "Amino Acid + Position", "Unimod"], how='outer') 

merged[['hPTM_ID', 'PTM Description', "Protein Accession", "Protein Description", "Position", "Amino Acid", "Amino Acid + Position", "Unimod","Normoxia vs 24hr Anoxia corrected p-value","Normoxia vs 4d Anoxia corrected p-value", "Normoxia vs Recovery corrected p-value", "4d Anoxia vs 24hr Anoxia corrected p-value", "Recovery vs 24hr Anoxia corrected p-value", "4d Anoxia vs Recovery corrected p-value"]]


merged.to_csv(base_path+'/SighPTMs.csv', index=False)
DataRecv24Filtered.to_csv(base_path+'/DataRecv24Filtered.csv', index=False)