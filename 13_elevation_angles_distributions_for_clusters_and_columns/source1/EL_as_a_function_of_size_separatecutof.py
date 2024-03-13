
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.stats import genextreme as gev


def plot_everything(datapath,ax,mycolor,legendlabel,cutoff):
	data=pd.read_csv(datapath+'spherical_coordinate.txt',sep='\t',header=None)
	#data=pd.read_csv('Embryo/spherical_coordinate12.txt',sep='\t',header=None)

	data=data.to_numpy()
	print(data.shape)

	factor=180/np.pi
	EL=data[:,2]*factor
	
	count=0
	total=0
	d={}

	clusize=[]
	el_fraction=[]
	
	cid=np.unique(data[:,0])
	print(len(cid))
	
	for k in range(len(cid)):
		index=np.where(data[:,0]==cid[k])
		#print(index[0])
		myEL=EL[index[0]]
		mydata=data[index[0]]
	
		temptotal=0
		tempcount=0
		dnoc={}
		for i in range(len(myEL)):
			name=str(mydata[i,4])+'#'+str(mydata[i,5])
			if name not in d:
				d[name]=1
				total+=1
				temptotal+=1
				if abs(myEL[i])>cutoff:
					count+=1
					tempcount+=1
					
			
			dnoc[mydata[i,4]]=1
			dnoc[mydata[i,5]]=1
		el_fraction.append(tempcount/temptotal)
		clusize.append(len(dnoc))	
							
				 		
						
	print('ful',total,count)	
	#print(clusize)
	#print(el_fraction)	
	
	X=np.unique(clusize)	
	Y=[]
	Z=[]
	Z0=[]
	for i in range(len(X)):
		temp=[]
		for j in range(len(clusize)):
			if clusize[j]==X[i]:
				temp.append(el_fraction[j])
		Y.append(np.mean(temp))
		Z.append(np.std(temp))
		Z0.append(0)
		
	#print(Y,Z)	
	stds=[Z0,Z]
		
	#fig, ax = plt.subplots()
	ax.bar(X, Y, yerr=stds, align='center', color='g', alpha=0.7, ecolor='black', error_kw=dict(lw=0.1, capsize=1, capthick=0.1),label=legendlabel)
	#ax.set_ylabel('Coefficient of Thermal Expansion ($\degree C^{-1}$)')
	#ax.set_xticks(x_pos)
	#ax.set_xticklabels(materials)
	#ax.set_title('Coefficent of Thermal Expansion (CTE) of Three Metals')
	#ax.yaxis.grid(True)

	# Save the figure and show
	
	return np.max(np.array(Z)+np.array(Y))
	
	
		

				
					
			

	












def main():
	datapath1='plothist_embryo/DF/'
	datapath2='plothist_postnatal/DF/'
	datapath3='plothist_embryo/PT/'
	datapath4='plothist_postnatal/PT/'

	cutoff=50

	fig,ax=plt.subplots(2,2,figsize=(8,5))
	

	a1=plot_everything(datapath1,ax[0,0],'r','E18.5 DF',cutoff)
	a2=plot_everything(datapath2,ax[0,1],'b','P40 DF',cutoff)
	a3=plot_everything(datapath3,ax[1,0],'r','E18.5 PT',cutoff)
	a4=plot_everything(datapath4,ax[1,1],'b','P40 PT',cutoff)
	ax[0,0].legend(prop={'size': 7},loc='upper right')
	ax[0,1].legend(prop={'size': 6},loc='upper right')
	ax[1,0].legend(prop={'size': 7},loc='upper right')
	ax[1,1].legend(prop={'size': 6},loc='upper right')
	ax[1,0].set_xlabel('cluster size (# of cells)')
	ax[1,1].set_xlabel('cluster size (# of cells)')
	ax[0,0].set_ylabel('% of doublets (|EL|> '+str(cutoff) +')')
	ax[1,0].set_ylabel('% of doublets (|EL|> '+str(cutoff) +')')
	largest=np.max([a1,a2,a3,a4])
	ax[0,0].set_ylim([0,largest])
	ax[0,1].set_ylim([0,largest])
	ax[1,0].set_ylim([0,largest])
	ax[1,1].set_ylim([0,largest])

	fig.tight_layout()
	fig.savefig('EL_as_a_function_of_size'+str(cutoff)+'.png',bbox_inches='tight',dpi=300)
	#fig.savefig('histogram_elevation_cluster_DF.png',bbox_inches='tight',dpi=300)
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
