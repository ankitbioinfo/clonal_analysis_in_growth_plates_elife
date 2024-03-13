

           

path=['S151_m2_distalfemur', 'S152_m3_distalfemur', 'S152_m4_distalfemur','S151_m2_proximaltibia', 'S152_m3_proximaltibia', 'S152_m4_proximaltibia']


def analysis(mypath):
	f=open(mypath+'/Cluster_NoDuplication_remove_nuclei.dat')
	#f=open(mypath+'/Cluster_NoDuplication.dat')
	cont=f.readlines()

	d={}
	d2={}
	for i in range(len(cont)):
		l=cont[i].split()
		for j in l:
			name=int(j)
			if name not in d:
				d[name]=1
				d2[name]=[i+1]
			else:
				d[name]+=1
				d2[name].append(i+1)
			
	s=[]			
	for key in d:
		if d[key]>1:
			tt=str(d2[key])
			print(key,d[key],d2[key])
			if tt not in s:
				s.append(tt)
				
	for i in range(len(s)):
		print(i,s[i])			
			
		
for i in range(6):
	mypath=path[i]
	print('\n\n',mypath)			
	analysis(mypath)
