#find the drug-drug interaction base on pubmed database
#convert generate to list
from Bio import Entrez 
import time
import random
import numpy as np
import scipy.stats
import math
import re
def getterms(myterm):
	start=myterm.find('(')
	return myterm[start+1:-1]
def all_lower(L1):
	return [s.lower() for s in L1]
def validate_drugname(drug):
	if re.match('^EC [0-9]+[0-9\.]+[0-9]+ \(.+\)$',drug) or re.match('^[A-Z0-9-]{9,11} \(.+\)$',drug):
		return True
	else:
		return False


def listsep(L1):
#in case the list being too long 
#sep 
	list_sep=list()
	i=0
	maxlen=9999	#max=9999
	len_list=len(L1)
	while (i+1)*maxlen<len_list:
		list_sep.append(L1[i*maxlen:(i+1)*maxlen])
		i+=1
		#print len_list
	list_sep.append(L1[i*maxlen:])
	return list_sep

def DDI_identify(L1):
# if the MH term (list) retated to ddi
	DDI_relatedMH=('drug interactions', 'drug agonism', 'drug partial agonism', 'drug antagonism',  'drug  inverse  agonism',  'drug  synergism',  'food-drug  interactions','herb-drug  interactions')
	for DDI_term in DDI_relatedMH:
		for  Mesh_term in all_lower(L1):
			if DDI_term in Mesh_term:
				return 1
	return 0


def article2drug(L1):
#convert the article based data(list) to drug based data(dict)
	drug_article=dict()
	for record in L1:
		if 'RN' in record:
			RN=record['RN']
			for drug in RN:
				if not drug.startswith('0'):
					if drug in drug_article:
						drug_article[drug]+=1
					else:
						drug_article[drug]=1
	return drug_article

def sepArticle(demo_drug):
#seperater the searched article based on the MeSh terms
#whether this article is ddi related
	#how many
	#print ">>>this is sepArticle"
	myterm=demo_drug
	Entrez.email = '812033546@qq.com'
	handle = Entrez.egquery(term=myterm)
	record = Entrez.read(handle)
	for row in record['eGQueryResult']:
		if row['DbName']=='pubmed':
			retmax=row['Count']
	print 'retmax',retmax#######################################
	#get idlist
	handle = Entrez.esearch(db='pubmed',term=myterm,retmax=retmax)##########
	record = Entrez.read(handle)
	idlist = record['IdList']
	#fetch the article
	from Bio import Medline
	idlist_sep=listsep(idlist)
	article_ddi=list();
	article_noneddi=list();
	for idlist0 in idlist_sep:
		handle = Entrez.efetch(db='pubmed', id=idlist0,rettype='medline',retmode='text')
		records=Medline.parse(handle) #records is a generater
		print '>>>convert generate to list', len(idlist0)
		records=list(records) #convert generate to list
		print '>>>done'
		for record in records:
			if 'MH' in record:#in case 
				MH=record['MH']
				if DDI_identify(MH):
					article_ddi.append(record)
				else:
					article_noneddi.append(record)
	return (article_ddi,article_noneddi)

def drug_pvalue(article_AB,demo_drug):
#caculate the p-value for the compouds
#the input is a tuple (len=2) contains two list of articles
	#print ">>>this is drug_pvalue"
	article_ddi=article_AB[0]
	article_noneddi=article_AB[1]
	count_articleDDI=len(article_ddi)
	count_articleNoneDDi=len(article_noneddi)
	print count_articleDDI,count_articleNoneDDi


	#calculate the p-value	
	drugs_ddi=article2drug(article_ddi)
	drugs_noneddi=article2drug(article_noneddi)
	drugs_pvalue=dict()
	enzymes_pvalue=dict()

	
	for drug in drugs_ddi:
		drug_frequency=drugs_ddi[drug]
		if drug in drugs_noneddi:
			drug_noneddi_frequency=drugs_noneddi[drug]
		else:
			drug_noneddi_frequency=0
		E_x=1.0*count_articleDDI*drug_noneddi_frequency/count_articleNoneDDi
		D_x=1.0*E_x*(1-1.0*drug_noneddi_frequency/count_articleNoneDDi)
		if drug_frequency>=4 and demo_drug not in  drug.lower() and validate_drugname(drug):
			if E_x!=0:
				pass
			else:
				E_x=1.0*count_articleDDI*0.1/count_articleNoneDDi
				D_x=E_x
			logp_value=scipy.stats.norm(E_x,D_x).logsf(drug_frequency)
			if logp_value<math.log(10**-6):
				if  not drug.startswith('EC') : 
					drugs_pvalue[drug]=(logp_value,drug_frequency,drug_noneddi_frequency)
				if drug.startswith('EC') :
					enzymes_pvalue[drug]=(logp_value,drug_frequency,drug_noneddi_frequency)
	return [drugs_pvalue,enzymes_pvalue,count_articleDDI,count_articleNoneDDi]

def main(demo_drug):
	print ">>>this is main"
	return drug_pvalue(sepArticle(demo_drug),demo_drug)
if __name__=="__main__":
	start=time.time()
	result=main('cytochrome p-450 enzyme system')
	drugs_pvalue=sorted(result[0].items(),key=lambda d:d[1])
	enzymes_pvalue=sorted(result[1].items(),key=lambda d:d[1])
	#print '>>>drugs\n' ,drugs_pvalue
	#print '>>>enzymes\n' ,enzymes_pvalue
	end=time.time()
	print ">>>Run time: %d s"%(end-start)


















