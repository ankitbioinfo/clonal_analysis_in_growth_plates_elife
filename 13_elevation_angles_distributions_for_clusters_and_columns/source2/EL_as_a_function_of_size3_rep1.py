
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.stats import genextreme as gev



def find_cluster_column_index(clusize,el_fraction):

	X=np.unique(clusize)	
	Y=[]
	Z=[]
	Z0=[]
	for i in range(len(X)):
		temp=[]
		for j in range(len(clusize)):
			if clusize[j]==X[i]:
				temp.append(el_fraction[j])
		a=np.mean(temp)
		b=np.mean(temp,axis=0)
		#print(a,b)		
		Y.append(np.nanmean(temp,axis=0))
		Z.append(np.std(temp,axis=0))
		Z0.append(0)
	stds=[Z0,Z]
	Y=np.array(Y)

	index1=[]
	index2=[]
	index3=[]
	index4=[]
	index5=[]
	index6=[]
	for i in range(len(X)):
		if X[i]<=5:
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
			
	a1=np.nanmean(Y[index1],axis=0)
	a2=np.nanmean(Y[index2],axis=0)
	a3=np.nanmean(Y[index3],axis=0)
	a4=np.nanmean(Y[index4],axis=0)
	a5=np.nanmean(Y[index5],axis=0)
	a6=np.nanmean(Y[index6],axis=0)
	YY=np.array([a1,a2,a3,a4,a5,a6])	
	return YY,Y,Z	


def plot_everything(datapath,PD_PC1_file,ax,mycolor,legendlabel,cutoff):

	PD_PC1={} 
	#0 means clusters and 90 means columns 
	
	for fi in range(len(PD_PC1_file)):
		fname=datapath+PD_PC1_file[fi]
		data=pd.read_csv(fname,sep=',',header=None)
		data=data.to_numpy()
		for i in range(data.shape[0]):
			#0 column is PD-PC1, 1 column is clone, 2nd column is section 
			name=str(fi+1)+'#'+str(i+1)  
			PD_PC1[name]=data[i,0]
		#print(data.shape,len(PD_PC1))
		

	data=pd.read_csv(datapath+'spherical_coordinate.txt',sep='\t',header=None)
	#data=pd.read_csv('Embryo/spherical_coordinate12.txt',sep='\t',header=None)

	data=data.to_numpy()
	#print(data.shape)

	factor=180/np.pi
	EL=data[:,2]*factor
	
	count=0
	total=0
	d={}

	clusize1=[]
	el_fraction1=[]
	clusize2=[]
	el_fraction2=[]
	
	newid=np.zeros((len(data)),dtype=object)
	for i in range(len(data)):
		newid[i]=str(int(data[i,0]))+'-'+str(data[i,4])[0]
	
	#cid=np.unique(data[:,0])
	cid=np.unique(newid)
	
	bininterval=[0,40,50,60,70,80,92]
	
	for k in range(len(cid)):
		index=np.where(newid==cid[k])
		#print(index[0])
		myEL=abs(EL[index[0]])
		mydata=data[index[0]]
	
		temptotal=0
		tempcount=np.zeros(len(bininterval)-1)
		dnoc={}
		for i in range(len(myEL)):
			name=str(mydata[i,4])+'#'+str(mydata[i,5])
			tname=str(mydata[i,4])[0]+'#'+str(int(mydata[i,0]))
			value=PD_PC1[tname]
			if name not in d:
					d[name]=1
					total+=1
					temptotal+=1
					for j in range(1,len(bininterval)):
						if (myEL[i]>=bininterval[j-1])&(myEL[i]<bininterval[j]) :
							count+=1
							tempcount[j-1]+=1
					
			dnoc[mydata[i,4]]=1
			dnoc[mydata[i,5]]=1


		if (value>=60)&(value<91):
			el_fraction1.append(tempcount/temptotal)
			clusize1.append(len(dnoc))	
		else:
			el_fraction2.append(tempcount/temptotal)
			clusize2.append(len(dnoc))		
							
				 		
						
	print('ful',total,count)	
	#print(clusize)
	print('ankur',len(el_fraction1),len(el_fraction2))	
	

	#print(len(index1)+len(index2)+len(index3)+len(index4)+len(index5)+len(index6),len(X))
	#print('1',Y.shape)
	#print('2',Y[index3].shape)
	
	YY,Y,Z=find_cluster_column_index(clusize1,el_fraction1)
	YY2,Y,Z=find_cluster_column_index(clusize2,el_fraction2)

	print('3',YY.shape,datapath,'column')
	
	for i in range(len(YY)):
		s=''
		for j in range(len(YY[0])):
			s+='\t'+str(YY[i][j])
		print(s)	
		
	print('3',YY.shape,datapath,'cluster')	
	for i in range(len(YY2)):
		s=''
		for j in range(len(YY2[0])):
			s+='\t'+str(YY2[i][j])
		print(s)	
	
	
	XX=['<=5','6-10','11-20','21-50','51-100','>100']
	
	
	ratio={'Col:0-40':YY[:,0],'Col:40:50':YY[:,1],'Col:50-60':YY[:,2],'Col:60-70':YY[:,3],'Col:70-80':YY[:,4],'Col:80-90':YY[:,5]}
	df1=pd.DataFrame(ratio,index=XX)
	my_colors=['y','c','g','m','r','b']
	df1.plot(ax=ax,kind='bar',stacked=True,width= 0.15,position=0.9,color=my_colors,alpha=1)

	
	ratio2={'Clu:0-40':YY2[:,0],'Clu:40:50':YY2[:,1],'Clu:50-60':YY2[:,2],'Clu:60-70':YY2[:,3],'Clu:70-80':YY2[:,4],'Clu:80-90':YY2[:,5]}
	df2=pd.DataFrame(ratio2,index=XX)
	df2.plot(ax=ax,kind='bar',stacked=True,width= 0.15, position=-0.5,color=my_colors,alpha=0.4)
	
	
	'''
	fig, ax = plt.subplots()
	df[['a', 'c']].plot.bar(stacked=True, width=0.1, position=1.5, colormap="bwr", ax=ax, alpha=0.7)
	df[['b', 'd']].plot.bar(stacked=True, width=0.1, position=-0.5, colormap="RdGy", ax=ax, alpha=0.7)
	df[['a', 'd']].plot.bar(stacked=True, width=0.1, position=0.5, colormap="BrBG", ax=ax, alpha=0.7)
	plt.legend(loc="upper center")
	plt.show()
	'''
	
	
	#ax.set_xticks(ind)
			
	#temp=np.floor(np.linspace(0,len(X),11))
	#temp=np.floor(np.linspace(0,len(X),11))

	ax.set_yticks([0,0.1,0.2,0.3,0.4,0.50,0.60,0.70,0.80,0.90,1.0])
	#ax.set_xticklabels([])
		
	#fig, ax = plt.subplots()
	#ax.bar(X, Y,  align='center', color='g', alpha=0.7, ecolor='black', error_kw=dict(lw=0.1, capsize=1, capthick=0.1),label=legendlabel) #yerr=stds,
	#ax.set_ylabel('Coefficient of Thermal Expansion ($\degree C^{-1}$)')
	#ax.set_xticks(x_pos)
	#ax.set_xticklabels(materials)
	ax.set_title(legendlabel)
	#ax.yaxis.grid(True)

	# Save the figure and show
	
	return np.max(np.array(Z)+np.array(Y))
	
	
		

				
					
			

	












def main():
	datapath1='plothist_embryo/DF/'
	datapath2='plothist_embryo_random/DF/rep1/'
	datapath3='plothist_embryo/PT/'
	datapath4='plothist_embryo_random/PT/rep1/'
	
	d1=['S153_m7_distalfemur.dat','S154_m3_distalfemur.dat','S154_m4_distalfemur.dat']
	d3=['S153_m7_proximaltibia.dat','S154_m3_proximaltibia.dat','S154_m4_proximaltibia.dat']

	d2=['S153_m7_distalfemur.dat','S154_m3_distalfemur.dat','S154_m4_distalfemur.dat']
	
	d4=['S153_m7_proximaltibia.dat', 'S154_m3_proximaltibia.dat', 'S154_m4_proximaltibia.dat']

	cutoff=0

	fig,ax=plt.subplots(2,2,figsize=(8,5))
	

	a1=plot_everything(datapath1,d1,ax[0,0],'r','E18.5 DF real',cutoff)
	a2=plot_everything(datapath2,d2,ax[0,1],'b','E18.5 DF random',cutoff)
	a3=plot_everything(datapath3,d3,ax[1,0],'r','E18.5 PT real',cutoff)
	a4=plot_everything(datapath4,d4,ax[1,1],'b','E18.5 PT random',cutoff)
	#ax[0,0].legend(prop={'size': 7},loc='upper right')
	ax[0,0].legend().set_visible(False)
	ax[0,1].legend(title="El angle",prop={'size': 6},loc='upper right',bbox_to_anchor=(1.28,1.02))
	#ax[1,0].legend(prop={'size': 7},loc='upper right')
	ax[1,0].legend().set_visible(False)

	ax[1,1].legend(title="El angle",prop={'size': 6},loc='upper right',bbox_to_anchor=(1.28,1.02))
	ax[1,0].set_xlabel('cluster size (# of cells)')
	ax[1,1].set_xlabel('cluster size (# of cells)')
	ax[0,0].set_ylabel('% of doublets')
	ax[1,0].set_ylabel('% of doublets')
	largest=np.max([a1,a2,a3,a4])
	
	#ax[1,0].legend(loc='lower left',fontsize=8, bbox_to_anchor=(0.7,1))
	
	
	ax[0,0].set_ylim([0.25,1.03])
	ax[0,1].set_ylim([0.25,1.03])
	ax[1,0].set_ylim([0.25,1.03])
	ax[1,1].set_ylim([0.25,1.03])
	

	fig.tight_layout()
	fig.savefig('EL_against_size_real_random_Rep1.png',bbox_inches='tight',dpi=300)
	#fig.savefig('histogram_elevation_cluster_DF.png',bbox_inches='tight',dpi=300)
	

	fig.savefig('EL_against_size_real_random_Rep1.svg',format='svg',bbox_inches='tight',dpi=1200)

	
	fig.clf()
	
main()


'''
import numpy as np
from matplotlib import pyplot

means   = [26.82,26.4,61.17,61.55]           # Mean Data 
stds    = [(0,0,0,0), [4.59,4.39,4.37,4.38]] # Standard deviation Data
peakval = ['26.82','26.4','61.17','61.55']   # String array of means

ind = np.arange(len(means))
width = 0.35
colours = ['red','blue','green','yellow']

pyplot.figure()
pyplot.title('Average Age')
pyplot.bar(ind, means, width, color=colours, align='center', yerr=stds, ecolor='k')
pyplot.ylabel('Age (years)')
pyplot.xticks(ind,('Young Male','Young Female','Elderly Male','Elderly Female'))

def autolabel(bars,peakval):
    for ii,bar in enumerate(bars):
        height = bars[ii]
        pyplot.text(ind[ii], height-5, '%s'% (peakval[ii]), ha='center', va='bottom')
#autolabel(means,peakval) 
pyplot.show()	
'''
