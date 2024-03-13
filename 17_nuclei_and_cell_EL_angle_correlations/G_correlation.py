import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.stats import genextreme as gev
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = 0.1 #set the value globally

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
plt.rc('font', family='Helvetica')

def plot_everything(datapath,ax,mycolor,legendlabel,alpha,lw):
    data=pd.read_csv(datapath+'spherical_coordinate.txt',sep='\t',header=None)
    #data=pd.read_csv('Embryo/spherical_coordinate12.txt',sep='\t',header=None)
    baddata=data.to_numpy()
    index=np.where(baddata[:,1]==1)
    data=baddata[index[0]]#only taking edges from contributing in the cluster
    clusterid=sorted(list(np.unique(data[:,0])))
    print(len(clusterid))


    clone_el={}
    for i in range(len(clusterid)):
        index=np.where(data[:,0]==clusterid[i])
        mydata=data[index[0],:]
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
            avgvalue+=d[key]/dval[key] #divide by because every nuclei comes two times

        factor=180/np.pi
        clone_el[clusterid[i]]=avgvalue*factor/len(d)
    return clone_el







def main():
	datapath1='cell/'
	datapath2='nuc/'


	fig,ax=plt.subplots(1,1,figsize=(3,1.8))


	cellel=plot_everything(datapath1,ax,'r','cell',0.3,5)
	nucel=plot_everything(datapath2,ax,'b','nuclei',0.6,1)

	x=[]
	y=[]
	for key in cellel:
		x.append(cellel[key])
		y.append(nucel[key])


	fig,ax=plt.subplots(1,1,figsize=(3,3))

	corr,_ = pearsonr(x,y)
	b, a = np.polyfit(x, y, deg=1)
	print(b,a)
	xseq = np.linspace(min(x), max(x), num=100)

	ax.plot(xseq, a + b * xseq, color="k", label='pearson r =%0.2f'%corr, lw=1);

	ax.plot(x,y,'bo',markerfacecolor="None")
	ax.legend( prop={'size': 6},loc='upper left')
	ax.set_title('n='+str(len(x)))
	ax.set_xlabel('Avg El of Clone [cell]')
	ax.set_ylabel('Avg EL of Clone [nuclei]')

	fig.savefig('EL_cell_nuclei_correlation.pdf',bbox_inches='tight',dpi=300)



main()
