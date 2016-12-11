import numpy as np 
import pandas as pd 
import re
#loda data from txt
input_dict_all=open('dict_all.txt','r')
dict_all=input_dict_all.read()
dict_all=eval(dict_all)
input_dict_all.close()

input_L_get=open('L_get.txt','r')
L_get=input_L_get.read()
L_get=eval(L_get)
input_L_get.close()

input_L_put=open('L_put.txt','r')
L_put=input_L_put.read()
L_put=eval(L_put)
input_L_put.close()

input_L_enzyme=open('L_enzyme.txt','r')
L_enzyme=input_L_enzyme.read()
L_enzyme=eval(L_enzyme)
input_L_enzyme.close()


def extract_drugName(drug_l):
	enzyme=re.findall('^EC [0-9]+[0-9\.]+[0-9]+ \((.+)\)$',drug_l)
	drug=re.findall('^[A-Z0-9-]{9,11} \((.+)\)$',drug_l)
	return (enzyme+drug)[0]

#creat a empty dataframe(1000*1000)
data_df=pd.DataFrame()

dict_drug=dict()
for drug in dict_all:
	dict_drug[drug]=dict_all[drug][0]
dict_drug1=dict()
for drug in dict_drug:
	dict_drug1[drug]=dict()
	for drug_l in dict_drug[drug]:
		drug_s=extract_drugName(drug_l)
		p_value=dict_drug[drug][drug_l][0]
		dict_drug1[drug][drug_s]=p_value
filecsv=open('result.csv','w')
for drug in dict_drug1:

	x=str(dict_drug1[drug].keys()).strip('[]\'')
	x=x.replace('\"','\'')
	x=x.strip('\'')
	y=drug+';'+x.replace('\', \'',';')+'\n'
	filecsv.writelines(y)
filecsv.close()



frame_drug=pd.DataFrame(dict_drug1)
frame_drug=frame_drug.fillna(0)
frame_drug.to_csv('ddi_matrix.csv')
