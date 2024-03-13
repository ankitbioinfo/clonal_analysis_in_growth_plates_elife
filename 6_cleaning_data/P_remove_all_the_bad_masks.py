

 

path=['S151_m2_distalfemur', 'S151_m2_proximaltibia','S152_m3_distalfemur', 'S152_m3_proximaltibia', 'S152_m4_distalfemur', 'S152_m4_proximaltibia']

remove_bad_masks=[[2,7,9,19],[80],[],[6,24,25,29,36,37,53],
[1,2,4,6,7,9,12,38,39,50,51,54,55,65,68,75,92,107,156,162,232,269,336,373,389,398],[634,1004,1271]]


def analysis(mypath,badmask):
	f=open(mypath+'/Cluster_NoDuplication_remove_nuclei.dat')
	cont=f.readlines()
	
	f=open(mypath+'/ClusterMask_NoDuplication.dat')
	mask=f.readlines()
	
	
	fw1=open(mypath+'/Cluster_good.dat','w')	
	fw2=open(mypath+'/ClusterMask_good.dat','w')
	
	for i in range(len(cont)):
		if (i+1) not in badmask:
			fw1.write(cont[i])
			fw2.write(mask[i]) 


for i in range(5,6):
	mypath=path[i]
	print('\n\n',mypath)			
	analysis(mypath,remove_bad_masks[i])	
