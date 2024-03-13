

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.stats import genextreme as gev


def find_cluster_column_index(clusize):
	X=np.unique(clusize)	
	Y=[]
	for i in range(len(X)):
		count=0
		for j in range(len(clusize)):
			if clusize[j]==X[i]:
				count=count+1 
		Y.append(count)
		
	Y=np.array(Y)	
	print('abc',len(clusize),sum(Y), len(X),len(Y))	
	
	index0=[]
	index1=[]
	index2=[]
	index3=[]
	index4=[]
	index5=[]
	index6=[]
	for i in range(len(X)):
		if X[i]<=2:
			index0.append(i)
		if (X[i]>2)&(X[i]<=5):
			index1.append(i)	
		if (X[i]>5)&(X[i]<=10):
			index2.append(i)
		if (X[i]>10)&(X[i]<=20):
			index3.append(i)
		if (X[i]>20)&(X[i]<=50):
			index4.append(i)
		if (X[i]>50)&(X[i]<=100):	
			index5.append(i)
		if X[i]>100:
			index6.append(i)
			
	a0=np.sum(Y[index0],axis=0)		
	a1=np.sum(Y[index1],axis=0)
	a2=np.sum(Y[index2],axis=0)
	a3=np.sum(Y[index3],axis=0)
	a4=np.sum(Y[index4],axis=0)
	a5=np.sum(Y[index5],axis=0)
	a6=np.sum(Y[index6],axis=0)
	YY=np.array([a0,a1,a2,a3,a4,a5,a6])
	
	return YY
							
		


def plot_everything(datapath,PD_PC1_file,ax,mycolor,legendlabel,cutoff):

	f=open(datapath+'degree_of_clusters.txt')
	clusize=[]
	for line in f:
		l=line.split(',')
		clusize.append(len(l)-1)	

	PD_PC1=[]
	#0 means clusters and 90 means columns 
	for fi in range(len(PD_PC1_file)):
		fname=datapath+PD_PC1_file[fi]
		data=pd.read_csv(fname,sep=',',header=None)
		data=data.to_numpy()
		for i in range(data.shape[0]):
			#0 column is PD-PC1, 1 column is clone, 2nd column is section 
			name=str(fi+1)+'#'+str(i+1)  
			PD_PC1.append(data[i,0])
		#print(data.shape,len(PD_PC1))
		

	print(len(clusize), len(PD_PC1))
	
	#print(clusize[0:10])
	
	clusize1=[]
	clusize2=[]
	for i in range(len(clusize)):
		value=PD_PC1[i]
		if (value>=60)&(value<91):
			clusize1.append(clusize[i])	
		else:
			clusize2.append(clusize[i])	
	

	a=find_cluster_column_index(clusize1)
	b=find_cluster_column_index(clusize2)
	
	#print('aaa', a, b)
	a=100*a/len(clusize)
	b=100*b/len(clusize)
	
	#print(a, b)
	


	myxlabel=['2','3','4','5','>5']
	myxlabel=['2','3-5','6-10','11-20','21-50','51-100','>100']

	ratio={'Cluster-like':b,'Column-like':a}
	df=pd.DataFrame(ratio,index=myxlabel)
	df.plot(ax=ax,kind='bar',stacked=True)
	ax.set_title(legendlabel)
	#ax[0][0].text(0,1,c[0])
	#ax[0][0].text(1,1,c[1])
	ax.set_ylabel('% frequency')
	#ax.legend(loc='lower left',fontsize=8, bbox_to_anchor=(0,1))
	ax.legend(loc='upper right',fontsize=8)













def main():
	datapath1='plothist_embryo/DF/'
	#datapath2='plothist_embryo_random/DF/rep1/'
	datapath2='plothist_embryo/PT/'
	#datapath4='plothist_embryo_random/PT/rep1/'
	
	
	d1=['1_S153_m7_distalfemur.dat','2_S154_m3_distalfemur.dat','3_S154_m4_distalfemur.dat']
	d2=['1_S153_m7_proximaltibia.dat','2_S154_m3_proximaltibia.dat','3_S154_m4_proximaltibia.dat']
	'''
	d2=['1_S151_m2_distalfemur.dat', '2_S152_m3_distalfemur.dat',  '3_S152_m4_distalfemur.dat']  
	#d4=['1_S151_m2_proximaltibia.dat','2_S152_m3_proximaltibia.dat', '3_S152_m4_proximaltibia.dat' ]
	d2=d1
	d4=d3
	'''
	
	#d2=['S153_m7_distalfemur.dat','S154_m3_distalfemur.dat','S154_m4_distalfemur.dat']
	#d4=['S153_m7_proximaltibia.dat', 'S154_m3_proximaltibia.dat', 'S154_m4_proximaltibia.dat']

	cutoff=0

	fig,ax=plt.subplots(1,2,figsize=(6,2.5))
	

	a1=plot_everything(datapath1,d1,ax[0],'r','E18.5 DF real',cutoff)
	a2=plot_everything(datapath2,d2,ax[1],'b','E18.5 PT real',cutoff)
	#a3=plot_everything(datapath3,d3,ax[1,0],'r','E18.5 PT real',cutoff)
	#a4=plot_everything(datapath4,d4,ax[1,1],'b','E18.5 PT random',cutoff)
	#ax[0,0].legend(prop={'size': 7},loc='upper right')
	#ax[0,0].legend().set_visible(False)
	#ax[0,1].legend(title="El angle",prop={'size': 6},loc='upper right',bbox_to_anchor=(1.39,1.03))
	#ax[1,0].legend(prop={'size': 7},loc='upper right')
	#ax[1,0].legend().set_visible(False)

	#ax[1,1].legend(title="El angle",prop={'size': 6},loc='upper right',bbox_to_anchor=(1.39,1.03))
	#ax[1,0].set_xlabel('cluster size (# of cells)')
	#ax[1,1].set_xlabel('cluster size (# of cells)')
	#ax[0,0].set_ylabel('% of doublets')
	#ax[1,0].set_ylabel('% of doublets')
	#largest=np.max([a1,a2,a3,a4])
	


	
	#ax[1,0].legend(loc='lower left',fontsize=8, bbox_to_anchor=(0.7,1))
	
	
	#ax[0,0].set_ylim([0.15,1.09])
	#ax[0,1].set_ylim([0.15,1.09])
	#ax[1,0].set_ylim([0.15,1.09])
	#ax[1,1].set_ylim([0.15,1.09])
	

	fig.tight_layout()
	fig.savefig('Freq_clones_real_imarismask.png',bbox_inches='tight',dpi=300)
	#fig.savefig('histogram_elevation_cluster_DF.png',bbox_inches='tight',dpi=300)
	fig.savefig('Freq_clones_real_imarismask.svg',format='svg',bbox_inches='tight',dpi=1200)
	fig.clf()
	
main()
