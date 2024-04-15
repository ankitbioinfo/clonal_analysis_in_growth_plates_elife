

filename=['/S153_m7_distalfemur','/S153_m7_proximaltibia','/S154_m3_distalfemur','/S154_m3_proximaltibia','/S154_m4_distalfemur','/S154_m4_proximaltibia']

id=5

f=open('micron15/merge_clones'+filename[id]+'/Cluster.dat')

clustermerged=[]
for line in f:
	a=[]
	l=line.split()
	for j in l:
		a.append(int(j))
	clustermerged.append(a)	
		
		
f=open('../MakeListNucleiLabelled'+filename[id]+'/Cluster_NoDuplication.dat')
clustermain=[]
for line in f:
	l=line.split()
	b=[]
	if len(l)>=2:
		for j in l:
			b.append(int(j))
	clustermain.append(b)	
	
	
fw=open('micron15/merge_clones'+filename[id]+'/ClusterMaskMerged.dat','w')	
for i in range(len(clustermerged)):
	node=clustermerged[i]
	found=[]
	for k in range(len(node)):
		for j in range(len(clustermain)):
			if node[k] in clustermain[j]:
				if (j+1) not in found:
					found.append(j+1)
					fw.write(str(j+1)+' ')
	#print(i+1,found)	
	fw.write('\n')			
				
				
		
		
	
		
		
print(len(a),len(b))
a=set(a)
b=set(b)
c=a.intersection(b)

print(len(c))	
print(b-a)		
