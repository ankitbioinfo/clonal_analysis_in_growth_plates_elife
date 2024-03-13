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

plt.rc('font', family='Helvetica')


def bimodal(x,mu1,mu2,s1,s2,p1):
    disA=np.exp(-0.5*(x-mu1)**2 / s1**2.) / (s1*np.sqrt(2.0 * np.pi))
    disB=np.exp(-0.5*(x-mu2)**2 / s2**2.) / (s2*np.sqrt(2.0 * np.pi))
    return p1*disA+(1-p1)*disB

def gumbel(x,mu,beta):
    return np.exp(-np.exp(-(x-mu)/beta))


def lognormal(x,mu,sigma):
	return (1.0/(x*sigma*np.sqrt(2*np.pi)))  *np.exp(- (np.log(x)-mu)**2 / (2*(sigma**2)) )


def exponential(x, a, b, c):
     return a*np.exp(-b*x) + c


def wigner(x,R):
	return (2/(np.pi*(R**2)))*np.sqrt(R**2 - x**2)


def plot_everything(datapath,ax,mycolor,legendlabel,alpha,lw):

	f=open(datapath+'Random_cluster.dat')

	deg=[]
	for line in f:
		l=line.split(',')
		for i in range(len(l)-1):
			deg.append(int(i))



	data=pd.read_csv(datapath+'spherical_coordinate.txt',sep='\t',header=None)
	#data=pd.read_csv('Embryo/spherical_coordinate12.txt',sep='\t',header=None)

	data=data.to_numpy()
	index=np.where(data[:,1]==1)
	data=data[index[0]]
	print(data.shape)

	factor=180/np.pi




	EL=data[:,3]*factor
	ax.hist(EL,bins=25,density=True,histtype='step',color=mycolor, facecolor=mycolor,alpha=1,label=legendlabel)
	ax.set_xlim([-95,95.1])
	ax.set_xlabel('Elevation')
	ax.set_ylabel('P(Elevation)')
	ax.set_xticks([-80,-40,0,40,80])
	#start, end = ax[1,0].get_xlim()
	#stepsize=30
	#ax[1,0].xaxis.set_ticks(np.arange(np.floor(start), end, stepsize))




	#ax[1,1].set_title('%5.3f exp(-%5.3fx) + %.2f' % tuple(popt))
	#ax[1,1].set_yscale('log')
	#ax[1,1].set_xscale('log')
	#ax[1,1].set_xlabel('deg')
	#ax[1,1].set_ylabel('P(deg)')
	#ax[1,1].set_xlim([0,40])

	#return legend2,lengend4


def main():
	datapath1='cell/'
	datapath2='nuc/'


	fig,ax=plt.subplots(1,1,figsize=(3,1.8))

	plot_everything(datapath1,ax,'r','cell',0.3,5)
	plot_everything(datapath2,ax,'b','nuclei',0.6,1)
	#ax[0,0].legend(prop={'size': 7},loc='lower center')
	#ax[0,1].legend(prop={'size': 6},loc='upper right')
	ax.legend(prop={'size': 5},loc='best')
	#ax[1,1].legend(prop={'size': 6},loc='upper right')

	#ax.set_xlim([-90,90])
	fig.tight_layout()
	#fig.savefig('Embryo_automated_cluster_cutoff_12.png',bbox_inches='tight',dpi=300)
	fig.savefig('histogram_cell_nuclei.pdf',bbox_inches='tight',dpi=300)
	#fig.savefig('histogram_cell_nuclei.pdf',format='svg',bbox_inches='tight',dpi=1200)
	fig.clf()

main()
