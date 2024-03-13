import pandas as pd 
import numpy as np 

df=pd.read_csv('matching_pairs_save.dat',header=None,sep='\t')
data=df.to_numpy()


fw=open('cell_nuc_pairs.dat','w')
n=5
cellid=sorted(list(np.unique(data[:,0])))
nucid=sorted(list(np.unique(data[:,1])))

cell2nuc={}
nuc2cell={}
alldata={}

for i in range(len(cellid)):
	index=np.where(data[:,0]==cellid[i])
	nuc=data[index[0],:]
	ind=np.argmin(nuc[:,2])
	neigh_nuc=nuc[ind,:]
	#print(i,neigh_nuc[1])
	cell2nuc[int(cellid[i])]=int(neigh_nuc[1]) 
	alldata[int(cellid[i])]=neigh_nuc 
	
for i in range(len(nucid)):
	index=np.where(data[:,1]==nucid[i])
	cel=data[index[0],:]
	ind=np.argmin(cel[:,2])
	neigh_cel=cel[ind,:]
	nuc2cell[int(nucid[i])]=int(neigh_cel[0])  	
	
good=0	
for kc in cell2nuc: 
	nid=cell2nuc[kc]
	#print(kc,nid)
	#for i in range(len(nid)):
	cid=nuc2cell[nid]
	#print(kc,nid,cid)	
	
	neigh_nuc =alldata[kc]
	pc1=neigh_nuc[3]
	#if pc1<45:
		#print(cellid[i],neigh_nuc)
	if kc==cid:	
		good+=1
		fw.write(str(kc)+'\t'+str(int(neigh_nuc[1]))+'\t'+str(neigh_nuc[2])+'\t'+str(neigh_nuc[3])+'\t'+str(neigh_nuc[4])+'\t'+str(neigh_nuc[5])+'\n')	
		
	
print(len(cell2nuc), good)		
		
		
		
