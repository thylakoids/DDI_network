import Queue
import medline3
import time
start=time.time()
start_now=start
start_write=start
dict_all=dict()
q_all=Queue.Queue()
L_enzyme=list()
L_put=list()
L_get=list()
#seed='cytochrome p-450 enzyme system'
seed='sertraline'

q_all.put(seed)
L_put.append(seed)
count=1
while not q_all.empty() and len(L_get)<1000:#############<1000
	myterm_name=q_all.get()
	print '>>>>>>>>>>>>>>>>>>>>>myterm:',count,myterm_name
	count=count+1
	L_get.append(myterm_name)

	###################################################### try except###############
	flag=1
	count_try=1
	while(flag):
		try:
			flag=0
			DDI_result=medline3.main(myterm_name)
		except Exception,e:
			flag=1
			count_try+=1
			print '+++++++++++++++++++++>>Warning:network error!'
			print e
			print 'Try time:',count_try,'(',count-1,myterm_name,')'
	#################################################################################
	dict_all[myterm_name]=DDI_result
	for tmp_d in list(DDI_result[0]):
		tmp_d=medline3.getterms(tmp_d)
		if tmp_d not in L_put:
			L_put.append(tmp_d)
			q_all.put(tmp_d)
	for tmp_e in list(DDI_result[1]):
		tmp_e=medline3.getterms(tmp_e)
		if tmp_e not in L_enzyme:
			L_enzyme.append(tmp_e)
	now=time.time()
	print '*********************time:',now-start_now,now-start_write,now-start
	if now-start_write>1800:
		output_dict_all=open('tmpdict_all.txt','w')
		output_dict_all.write(str(dict_all))
		output_dict_all.close()
		
		output_L_get=open('tmpL_get.txt','w')
		output_L_get.write(str(L_get))
		output_L_get.close()
		
		output_L_put=open('tmpL_put.txt','w')
		output_L_put.write(str(L_put))
		output_L_put.close()
		
		output_L_enzyme=open('tmpL_enzyme.txt','w')
		output_L_enzyme.write(str(L_enzyme))
		output_L_enzyme.close()
		start_write=now
	start_now=now

output_dict_all=open('dict_all.txt','w')
output_dict_all.write(str(dict_all))
output_dict_all.close()

output_L_get=open('L_get.txt','w')
output_L_get.write(str(L_get))
output_L_get.close()

output_L_put=open('L_put.txt','w')
output_L_put.write(str(L_put))
output_L_put.close()

output_L_enzyme=open('L_enzyme.txt','w')
output_L_enzyme.write(str(L_enzyme))
output_L_enzyme.close()



	
