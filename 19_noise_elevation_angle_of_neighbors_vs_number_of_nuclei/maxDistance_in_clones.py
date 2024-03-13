
import pandas as pd 
import numpy as np 
import networkx as nx 
import matplotlib.pyplot as plt


def el_angle_of_nodes(mydata):
        factor=180/np.pi
        d={}
        dval={}
        for j in range(len(mydata)):
            d[mydata[j,5]]=0
            d[mydata[j,6]]=0
            dval[mydata[j,5]]=0
            dval[mydata[j,6]]=0

        for j in range(len(mydata)):
            #print(mydata[j,3])
            d[mydata[j,5]]=d[mydata[j,5]]+abs(mydata[j,3])
            d[mydata[j,6]]=d[mydata[j,6]]+abs(mydata[j,3])
            dval[mydata[j,5]]=dval[mydata[j,5]]+1
            dval[mydata[j,6]]=dval[mydata[j,6]]+1

        avgvalue=0
        for key in d:
            #print(key,d[key],dval[key])
            #avgvalue+=d[key]/dval[key] #divide by because every nuclei comes two times
            d[key]=factor*d[key]/dval[key] #divide by because every nuclei comes two times


        #clone_el[clusterid[i]]=avgvalue*factor/len(d)
        return d 


def main_analysis(embryo,myfile):
	
	clonesize={}

	for fi in range(len(myfile)):
		
	
	
		df=pd.read_excel(embryo+myfile[fi]+'/nuclei_column_stats_data.xlsx')
		data=df['AngleBetweenClusterPC1AndBone_PD']
		pdpc1=data.to_numpy()
		#print(pdpc1)
	
		df=pd.read_csv(embryo+myfile[fi]+'/spherical_coordinate.txt',sep='\t',header=None)
		data=df.to_numpy()
		el_participate=np.where(data[:,1]==1)
		data=data[el_participate[0],:]
		cloneid=np.unique(data[:,0])
		no_of_clone=len(cloneid)
		
		totalnuc=0
		for i in range(len(cloneid)):
			index=np.where(data[:,0]==cloneid[i])
			clone_data=data[index[0],:]
			avg_el=el_angle_of_nodes(clone_data)
			edges=clone_data[:,[5,6]]
			dist=2*clone_data[:,4]
			G=nx.Graph()
			G.add_edges_from(edges)
			value=G.degree()
			deg1 = [d for n, d in G.degree()]  # degree sequence
			deg=[]
			EL=[]
			totalnuc+=len(deg1)
			#print(i,deg1,max(dist))
			clonesize[len(deg1)]=max(dist)
			
		print(len(cloneid),len(clonesize))	
		
	mu=[]
	std=[]
	deg=[0,5,10,20,50,100,1000]
	
	spread=[]
	for j in range(1,len(deg)):
		spread.append([])
	
	for key in clonesize:
		for j in range(1,len(deg)):
			if deg[j-1]<key<=deg[j]:
				spread[j-1].append(clonesize[key])
				
	
	for i in range(len(spread)):
		mu.append(np.mean(spread[i]))
		std.append(np.std(spread[i]))			
					
		
	return mu,std	
	

def plot_things(savename,mu,std,titlename):
	x=[1,2,3,4,5,6]
	mu=np.array(mu)
	std=np.array(std)
	fig,ax=plt.subplots(1,1,figsize=(4,3))
	ax.plot(x,mu,'bo-',ms=1.5,label='<EL> of nuclei in Columns')	
	#ax.plot(x,std,'ro',ms=0.5,label='<EL> of nuclei in Cluster')
	
	
	ax.fill_between(x, mu-std, mu+std,color="green", alpha=0.25)#, interpolate=True,  label="Positive")

	ax.set_xticks(x)
	xlabel=['2<=size<=5','5<size<=10','10<size<=20','20<size<=50','50<size<=100','size>100']  
	ax.set_xticklabels(xlabel,rotation=30)

	ax.set_xlabel('Clone sizes',fontsize=9)
	ax.set_ylabel('farthest apart doublet[micron]',fontsize=9)
	ax.set_title(titlename)
	#ax.set_xlim([0,1+value])
		

	fig.tight_layout()
	#fig.savefig('Embryo_automated_cluster_cutoff_12.png',bbox_inches='tight',dpi=300)
	#fig.savefig(savename,bbox_inches='tight',dpi=300)
	fig.savefig(savename+'.svg',format='svg',bbox_inches='tight',dpi=1200)



def main():
	df_file=['S153_m7_distalfemur','S154_m3_distalfemur','S154_m4_distalfemur']
	pt_file=['S153_m7_proximaltibia','S154_m3_proximaltibia', 'S154_m4_proximaltibia']
	embryo='RealSampling_embryo/'

	mu,std=main_analysis(embryo,df_file)
	plot_things('max_distance_embryo_DF',mu,std,'Embryo DF')

	
	mu,std=main_analysis(embryo,pt_file)
	plot_things('max_distance_embryo_PT',mu,std,'Embryo PT')
	
	P40='PT_MakeListNucleiLabelled/'
	df_file=['S151_m2_distalfemur','S152_m3_distalfemur','S152_m4_distalfemur']
	pt_file=['S151_m2_proximaltibia','S152_m3_proximaltibia','S152_m4_proximaltibia']
	
	
	mu,std=main_analysis(P40,df_file)
	plot_things('max_distance_P40_DF',mu,std,'Postnatal DF')

	mu,std=main_analysis(P40,pt_file)
	plot_things('max_distance_P40_PT',mu,std,'Postnatal PT')
	


main()


