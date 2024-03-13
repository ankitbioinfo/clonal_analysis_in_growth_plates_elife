
import matplotlib.pyplot as plt 
import numpy as np


def plot_figure(pd_pc1_df,ax,mycolor,legendlabel,subplot,design,alpha):
	print(len(pd_pc1_df))
	ax[subplot].hist(pd_pc1_df,bins=30,density=True,histtype=design, color=mycolor, facecolor=mycolor,alpha=alpha,label=legendlabel)

def readfiles(path,df_file):
	pd_pc1_df=[]
	for fi in range(len(df_file)):
		f=open(path+df_file[fi],'r')
		cont=f.readlines()
		for i in range(len(cont)):
			l=cont[i].split(',')
			pd_pc1_df.append(float(l[0]))
			
	return pd_pc1_df 
	

def main():
	path2='15_micron_2colorCoded3d_for_PC1_and_PD_merged_mask/'
	path1='2colorCoded3d_for_PC1_and_PD_embryo_uniclone/'
	path3='2colorCoded3d_for_PC1_and_PD_postnatal/'

				
	fig,ax=plt.subplots(1,2,figsize=(8,3))
	
	df_file=['S153_m7_distalfemur.dat', 'S154_m3_distalfemur.dat',  'S154_m4_distalfemur.dat']
	pt_file=['S153_m7_proximaltibia.dat', 'S154_m3_proximaltibia.dat',  'S154_m4_proximaltibia.dat']

	pd_pc1_df=readfiles(path1,df_file)
	pd_pc1_pt=readfiles(path1,pt_file)
	plot_figure(pd_pc1_df,ax,'b','E18.5 DF uniclone',0,'step',1)		
	plot_figure(pd_pc1_pt,ax,'b','E18.5 PT uniclone',1,'step',1)


	pd_pc1_df=readfiles(path2,df_file)
	pd_pc1_pt=readfiles(path2,pt_file)
	plot_figure(pd_pc1_df,ax,'g','E18.5 DF merge clone',0,'step',1)		
	plot_figure(pd_pc1_pt,ax,'g','E18.5 PT merge clone',1,'step',1)





	df_file=['S151_m2_distalfemur.dat','S152_m3_distalfemur.dat','S152_m4_distalfemur.dat']
	pt_file=['S151_m2_proximaltibia.dat','S152_m3_proximaltibia.dat','S152_m4_proximaltibia.dat']

	pd_pc1_df=readfiles(path3,df_file)
	pd_pc1_pt=readfiles(path3,pt_file)
	plot_figure(pd_pc1_df,ax,'r','P40 DF uniclone',0,'stepfilled',0.2)		
	plot_figure(pd_pc1_pt,ax,'r','P40 PT uniclone',1,'stepfilled',0.2)


	#ax[0,0].set_xlim([-180,180.1])
	ax[0].set_xlabel('[Angle between PD and PC1]')
	ax[1].set_xlabel('[Angle between PD and PC1]')
	ax[0].set_ylabel('P(PD-PC1)')
	ax[0].legend(prop={'size': 6},loc='upper right')
	ax[1].legend(prop={'size': 6},loc='upper right')
	#	start, end = ax[0,0].get_xlim()
	#	stepsize=60
	ax[0].xaxis.set_ticks(np.arange(0, 90.1, 10))
	ax[1].xaxis.set_ticks(np.arange(0, 90.1, 10))
	ax[1].set_ylim([0,0.09])
	ax[0].set_ylim([0,0.09])
				
	fig.tight_layout()
	#fig.savefig('Embryo_automated_cluster_cutoff_12.png',bbox_inches='tight',dpi=300)
	fig.savefig('histogram_PD_PC1.png',bbox_inches='tight',dpi=300)
	fig.clf()
				
				
main()				
