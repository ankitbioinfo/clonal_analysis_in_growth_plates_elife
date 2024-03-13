
import pandas as pd 
import numpy as np 
import networkx as nx 
import matplotlib.pyplot as plt


def el_angle_of_nodes(mydata):
        factor=180/np.pi
        d={}
        dval={}
        for j in range(len(mydata)):
            d[mydata[j,5]]=[]
            d[mydata[j,6]]=[]
            dval[mydata[j,5]]=0
            dval[mydata[j,6]]=0

        for j in range(0,len(mydata),2):
            #print(mydata[j,:])
            d[mydata[j,5]].append(factor*abs(mydata[j,3]))
            d[mydata[j,6]].append(factor*abs(mydata[j,3]))
            dval[mydata[j,5]]=dval[mydata[j,5]]+1
            dval[mydata[j,6]]=dval[mydata[j,6]]+1

        #avgvalue=0
        #for key in d:
            #print(key,d[key],dval[key])
            #avgvalue+=d[key]/dval[key] #divide by because every nuclei comes two times
            #d[key]=factor*d[key]/dval[key] #divide by because every nuclei comes two times


        #clone_el[clusterid[i]]=avgvalue*factor/len(d)
        return d 


def main_analysis(embryo,myfile):
	xcol=[]
	ycol=[]
	xclu=[]
	yclu=[]
	
	savedegcol={}
	savedegclu={}
	for i in range(1,500):
		savedegcol[i]=[0,0]
		savedegclu[i]=[0,0]

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
		
		anglecount=[0,60,91]
		#anglecount=[0,15,30,45,60,75,91]
		
		
		for i in range(len(cloneid)):
			index=np.where(data[:,0]==cloneid[i])
			clone_data=data[index[0],:]
			avg_el=el_angle_of_nodes(clone_data)
			edges=clone_data[:,[5,6]]
			G=nx.Graph()
			G.add_edges_from(edges)
			value=G.degree()
			deg1 = [d for n, d in G.degree()]  # degree sequence
			deg=[]
			EL=[]
			#print('\n\n',i,len(value),value,avg_el)
			print(avg_el)
			
			
		print(no_of_clone)
	return savedegcol, savedegclu	
	

def plot_things(savename,xcol,xclu,every_nth):
	fig,ax=plt.subplots(2,1,figsize=(5,4))
	my_colors=['y','c','g','m','r','b']
	mylabel=['0-<15','15-<30','30-<45','45-<60','60-<75','75-90']
	header=[xcol,xclu]
	
	for k in range(len(header)):
		ratio=[]
		myratio={}
		XX=[]
		for i in range(1,50):
			p=sum(header[k][i])
			l=np.array(header[k][i])/p
			if p!=0:
				ratio.append(l)
				XX.append(i)
		ratio=np.array(ratio)
		#print(ratio.shape)

		for i in range(ratio.shape[1]):
		 	myratio[mylabel[i]]=ratio[:,i]
		
		
		df1=pd.DataFrame(myratio,index=XX)
		df1.plot(ax=ax[k],kind='bar',stacked=True,width= 0.25,color=my_colors,alpha=1,rot=0)
		ax[k].set_xlabel('number of nuclei neighbors')
		ax[k].set_ylabel('El angle of neighbors')
		#ax[k].tick_params(axis='x', labelsize=8)
		
		 
		for n, label in enumerate(ax[k].xaxis.get_ticklabels()):
			if n % every_nth != 0:
				label.set_visible(False)

	
	#ax.bar(xcol,ycol,'bs',ms=1.5,label='Columns')	
	#ax.plot(xclu,yclu,'ro',ms=0.5,label='Cluster')	
		
	#ax.legend(prop={'size': 7},loc='upper right')
	
	ax[0].legend(title="Columns",prop={'size': 6},loc='upper right',bbox_to_anchor=(1.25,1.03))
	ax[1].legend(title="Clusters",prop={'size': 6},loc='upper right',bbox_to_anchor=(1.25,1.03))


		

	fig.tight_layout()
	#fig.savefig('Embryo_automated_cluster_cutoff_12.png',bbox_inches='tight',dpi=300)
	#fig.savefig(savename,bbox_inches='tight',dpi=300)
	fig.savefig(savename+'.svg',format='svg',bbox_inches='tight',dpi=1200)



def main():
	df_file=['S153_m7_distalfemur','S154_m3_distalfemur','S154_m4_distalfemur']
	pt_file=['S153_m7_proximaltibia','S154_m3_proximaltibia', 'S154_m4_proximaltibia']
	embryo='RealSampling_embryo/'

	xcol,xclu=main_analysis(embryo,df_file)
	#plot_things('nuclei_neighbors_embryo_DF',xcol,xclu,2)

	
	xcol,xclu=main_analysis(embryo,pt_file)
	#plot_things('nuclei_neighbors_embryo_PT',xcol,xclu,2)
	
	P40='PT_MakeListNucleiLabelled/'
	df_file=['S151_m2_distalfemur','S152_m3_distalfemur','S152_m4_distalfemur']
	pt_file=['S151_m2_proximaltibia','S152_m3_proximaltibia','S152_m4_proximaltibia']
	
	
	xcol,xclu=main_analysis(P40,df_file)
	#plot_things('nuclei_neighbors_P40_DF',xcol,xclu,4)

	xcol,xclu=main_analysis(P40,pt_file)
	#plot_things('nuclei_neighbors_P40_PT',xcol,xclu,4)
	


main()


