
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
	print('\n\n')
	xcol=[]
	ycol=[]
	xclu=[]
	yclu=[]
	savedegcol={}
	savedegclu={}
	for i in range(1,500):
		savedegcol[i]=0
		savedegclu[i]=0	

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
			G=nx.Graph()
			G.add_edges_from(edges)
			value=G.degree()
			deg1 = [d for n, d in G.degree()]  # degree sequence
			deg=[]
			EL=[]
			totalnuc+=len(deg1)
			for n,d in G.degree():
				deg.append(d)
				#print(n,d,avg_el)
				EL.append(avg_el[n])

			
			#print(i+1, len(deg), deg1,deg,EL)
			
			for k in range(len(deg1)):
				if pdpc1[i]>60:
					xcol.append(deg1[k])
					ycol.append(EL[k])
					savedegcol[deg1[k]]+=1
				else:
					xclu.append(deg1[k])
					yclu.append(EL[k])
					savedegclu[deg1[k]]+=1	
		print(no_of_clone,'total',totalnuc,sum(savedegclu.values()), sum(savedegcol.values()))
			
		
		
	return xcol,ycol,savedegcol,	xclu,yclu,savedegclu	
	
def find_prop(tcol,pcol,cutoff):

	totnuc_col=sum(pcol.values())
	print('t',totnuc_col,cutoff)

	prop_col=[]	
	lessdn5_col=0
	for i in range(len(tcol)):
		prop_col.append(100*pcol[tcol[i]]/totnuc_col)
		if tcol[i]<=cutoff:
			lessdn5_col+=prop_col[i]
	return prop_col,lessdn5_col		
		
	
	

def plot_things(savename,xcol,ycol,pcol,xclu,yclu,pclu,age):
	fig,ax=plt.subplots(1,1,figsize=(5,3))
	ax.plot(xcol,ycol,'bs',ms=1.5,label='<EL> of nuclei in Columns')	
	ax.plot(xclu,yclu,'ro',ms=0.5,label='<EL> of nuclei in Cluster')
	
	
	#print(totnuc_col, totnuc_clu)
	tcol=sorted(list(np.unique(xcol)))
	tclu=sorted(list(np.unique(xclu)))


	prop_col,lessdn5_col=find_prop(tcol,pcol,5)
	prop_col,lessdn10_col=find_prop(tcol,pcol,10)

	
	prop_clu,lessdn5_clu=find_prop(tclu,pclu,5)
	prop_clu,lessdn10_clu=find_prop(tclu,pclu,10)
	
	if age==1:
		prop_col,lessdn20_col=find_prop(tcol,pcol,20)
		prop_clu,lessdn20_clu=find_prop(tclu,pclu,20)
		print('fraction', lessdn5_col,lessdn10_col,lessdn20_col,'\t',lessdn5_clu,lessdn10_clu, lessdn20_clu)	
	
		
		
	value=	max([max(tcol),max(tclu)])
		
	print('fraction', lessdn5_col,lessdn10_col,'\t',lessdn5_clu,lessdn10_clu)	
		
	ax.plot([5.5,5.5], [0,90],'k:',lw=0.5)	
	ax.text(2.5,91,'%0.1f'%lessdn5_col+'%',color='b',fontsize=5)	
	ax.text(2.5,88,'%0.1f'%lessdn5_clu+'%',color='r',fontsize=5)
	
	
	if age==0:
		ax.plot([10.5,10.5], [0,90],'k:',lw=0.5)	
		ax.text(7.5,91,'%0.1f'%(lessdn10_col-lessdn5_col)+'%',color='b',fontsize=5)	
		ax.text(7.5,88,'%0.1f'%(lessdn10_clu-lessdn5_clu)+'%',color='r',fontsize=5)	
	
	
	if age==1:
		ax.plot([10.5,10.5], [0,90],'k:',lw=0.5)	
		ax.text(10.5,91,'%0.1f'%(lessdn10_col-lessdn5_col)+'%',color='b',fontsize=5)	
		ax.text(10.5,88,'%0.1f'%(lessdn10_clu-lessdn5_clu)+'%',color='r',fontsize=5)	
		
		ax.plot([20.5,20.5], [0,90],'k:',lw=0.5)	
		ax.text(20.5,91,'%0.1f'%(lessdn20_col-lessdn10_col)+'%',color='b',fontsize=5)	
		ax.text(20.5,88,'%0.1f'%(lessdn20_clu-lessdn10_clu)+'%',color='r',fontsize=5)	
	
		
	
	#ax.plot([0,1+value], [5,5],'k:',lw=0.5)	
	#ax.plot([0,1+value], [10,10],'k:',lw=0.5)
	#ax.plot([0,1+value], [1,1],'k:',lw=0.5)		
		
	#ax.plot(tcol,prop_col,'m-',lw=1,label='nuclei prop only in Column')	
	#ax.plot(tclu,prop_clu,'g--',lw=1,label='nuclei prop only in Cluster')	
		
	ax.legend(prop={'size': 7},loc='upper right')

	ax.set_xlabel('nuclei neighbors')
	ax.set_ylabel('avg El angle of nuclei')
	ax.set_xlim([0,1+value])
	ax.set_ylim([0,95])
		

	fig.tight_layout()
	#fig.savefig('Embryo_automated_cluster_cutoff_12.png',bbox_inches='tight',dpi=300)
	#fig.savefig(savename,bbox_inches='tight',dpi=300)
	fig.savefig(savename+'.svg',format='svg',bbox_inches='tight',dpi=1200)



def main():
	df_file=['S153_m7_distalfemur','S154_m3_distalfemur','S154_m4_distalfemur']
	pt_file=['S153_m7_proximaltibia','S154_m3_proximaltibia', 'S154_m4_proximaltibia']
	embryo='RealSampling_embryo/'

	'''
	xcol,ycol,pcol,xclu,yclu,pclu=main_analysis(embryo,df_file)
	plot_things('deg_vs_el_embryo_DF',xcol,ycol,pcol,xclu,yclu,pclu,0)

	xcol,ycol,pcol,xclu,yclu,pclu=main_analysis(embryo,pt_file)
	plot_things('deg_vs_el_embryo_PT',xcol,ycol,pcol,xclu,yclu,pclu,0)
	
	'''
	P40='PT_MakeListNucleiLabelled/'
	df_file=['S151_m2_distalfemur','S152_m3_distalfemur','S152_m4_distalfemur']
	pt_file=['S151_m2_proximaltibia','S152_m3_proximaltibia','S152_m4_proximaltibia']
	
	
	xcol,ycol,pcol,xclu,yclu,pclu=main_analysis(P40,df_file)
	plot_things('deg_vs_el_P40_DF',xcol,ycol,pcol,xclu,yclu,pclu,1)

	xcol,ycol,pcol,xclu,yclu,pclu=main_analysis(P40,pt_file)
	plot_things('deg_vs_el_P40_PT',xcol,ycol,pcol,xclu,yclu,pclu,1)
	


main()


