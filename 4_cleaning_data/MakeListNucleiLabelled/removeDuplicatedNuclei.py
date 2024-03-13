
            

path=['S151_m2_distalfemur', 'S152_m3_distalfemur', 'S152_m4_distalfemur','S151_m2_proximaltibia', 'S152_m3_proximaltibia', 'S152_m4_proximaltibia']

mypath=path[5]

def analysis(mypath):

	f=open(mypath+'/Cluster_NoDuplication.dat')
	cont=f.readlines()

	d={}
	d2={}
	data=[]
	for i in range(len(cont)):
		l=cont[i][0:-1].split()
		data.append(set(l))
		for j in l:
			name=int(j)
			if name not in d:
				d[name]=1
				d2[name]=[i]
			else:
				d[name]+=1
				d2[name].append(i)
			
	s=[]			
	for key in d:
		if d[key]>1:
			n=len(d2[key])
			if n==2:
				t1=min(d2[key])
				t2=max(d2[key])
				p=data[t1].difference(data[t2])
				data[t1]=p 
			else:
				t=sorted(d2[key])
				#print('a1',len(data[t1]))
				#print('a2',len(data[t2]))
				for j in range(1,len(t)):
					p=data[t[0]].difference(data[t[j]])
					data[t[0]]=p 
					
				for j in range(2,len(t)):
					p=data[t[1]].difference(data[t[j]])
					data[t[1]]=p 	
					
				for j in range(3,len(t)):
					p=data[t[2]].difference(data[t[j]])
					data[t[2]]=p 	
					
				for j in range(4,len(t)):
					p=data[t[3]].difference(data[t[j]])
					data[t[3]]=p 
					
				for j in range(5,len(t)):
					p=data[t[4]].difference(data[t[j]])
					data[t[4]]=p 				
					
			
		
	f=open(mypath+'/Cluster_NoDuplication_remove_nuclei.dat','w')		
	for i in range(len(data)):
		l=sorted(list(data[i]))
		for j in range(len(l)):
			f.write(str(l[j])+' ')
		f.write('\n')	
			
			
	
for i in range(6):
	mypath=path[i]
	print('\n\n',mypath)			
	analysis(mypath)			
				
