#Author: Chelsea Hughes

import numpy as np
from pandas import read_csv
import pandas as pd
from numpy import inf

base_path = "/Users/chelseahughes/Desktop/Histone Analysis/Calculations for WS40NE Samples/"

#Making Venn diagram for all sig PTMs
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
from venn import venn

sigdata=read_csv(base_path+'/SighPTMs.csv',header=0)
sigdata=sigdata[sigdata.columns[8:14]]
for index, row in sigdata.iterrows():
    for col in sigdata:
        if not pd.isna(row[col]):
            row[col] = index
dictionaryInstance = sigdata.to_dict(orient="list")
cleaned_dict = {}
for data in dictionaryInstance:
    nan_cleared_list = [x for x in dictionaryInstance[data] if x==x]
    cleaned_dict[data]=set(nan_cleared_list)
for data in dictionaryInstance:
    dictionaryInstance[data]=set(dictionaryInstance[data])
x = venn(cleaned_dict, cmap="plasma", figsize=(40,40), fontsize=30, legend_loc="upper left")
#plt.legend(x.figure[0], bbox_to_anchor=(1,0.5), loc="center right", fontsize=10, 
           #bbox_transform=plt.gcf().transFigure)
#plt.subplots_adjust(left=0.0, bottom=0.1, right=0.45)
#x.figure.legend = x.figure.legend(loc="center right")
x.figure.savefig("TotalSighPTMs.jpeg", bbox_inches='tight',pad_inches=1, dpi=1000)
#venn2([set(lst1), set(lst2)], set_labels = ('Drug A responsive genes', 'Drug B responsive genes'))
#plt.title('Comparison of Significant hPTMs\n')

sigdata=read_csv(base_path+'/SighPTMsH1.csv',header=0)
sigdata=sigdata[sigdata.columns[8:14]]
for index, row in sigdata.iterrows():
    for col in sigdata:
        if not pd.isna(row[col]):
            row[col] = index
dictionaryInstance = sigdata.to_dict(orient="list")
cleaned_dict = {}
for data in dictionaryInstance:
    nan_cleared_list = [x for x in dictionaryInstance[data] if x==x]
    cleaned_dict[data]=set(nan_cleared_list)
x = venn(cleaned_dict, cmap="plasma", figsize=(40,40), fontsize=30, legend_loc="upper left")
x.figure.savefig("SighPTMsH1.jpeg", bbox_inches='tight',pad_inches=1, dpi=1000)

sigdata=read_csv(base_path+'/SighPTMsH3.csv',header=0)
sigdata=sigdata[sigdata.columns[8:14]]
for index, row in sigdata.iterrows():
    for col in sigdata:
        if not pd.isna(row[col]):
            row[col] = index
dictionaryInstance = sigdata.to_dict(orient="list")
cleaned_dict = {}
for data in dictionaryInstance:
    nan_cleared_list = [x for x in dictionaryInstance[data] if x==x]
    cleaned_dict[data]=set(nan_cleared_list)
x = venn(cleaned_dict, cmap="plasma", figsize=(40,40), fontsize=30, legend_loc="upper left")
x.figure.savefig("SighPTMsH3.jpeg", bbox_inches='tight',pad_inches=1, dpi=1000)

sigdata=read_csv(base_path+'/SighPTMsH2A.csv',header=0)
sigdata=sigdata[sigdata.columns[8:14]]
for index, row in sigdata.iterrows():
    for col in sigdata:
        if not pd.isna(row[col]):
            row[col] = index
dictionaryInstance = sigdata.to_dict(orient="list")
cleaned_dict = {}
for data in dictionaryInstance:
    nan_cleared_list = [x for x in dictionaryInstance[data] if x==x]
    cleaned_dict[data]=set(nan_cleared_list)
x = venn(cleaned_dict, cmap="plasma", figsize=(40,40), fontsize=30, legend_loc="upper left")
x.figure.savefig("SighPTMsH2A.jpeg", bbox_inches='tight',pad_inches=1, dpi=1000)

sigdata=read_csv(base_path+'/SighPTMsH2B.csv',header=0)
sigdata=sigdata[sigdata.columns[8:14]]
for index, row in sigdata.iterrows():
    for col in sigdata:
        if not pd.isna(row[col]):
            row[col] = index
dictionaryInstance = sigdata.to_dict(orient="list")
cleaned_dict = {}
for data in dictionaryInstance:
    nan_cleared_list = [x for x in dictionaryInstance[data] if x==x]
    cleaned_dict[data]=set(nan_cleared_list)
x = venn(cleaned_dict, cmap="plasma", figsize=(40,40), fontsize=30, legend_loc="upper left")
x.figure.savefig("SighPTMsH2B.jpeg", bbox_inches='tight',pad_inches=1, dpi=1000)

