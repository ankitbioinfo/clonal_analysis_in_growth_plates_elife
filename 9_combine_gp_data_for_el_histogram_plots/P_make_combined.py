      


df_file=['S151_m2_distalfemur','S152_m3_distalfemur','S152_m4_distalfemur']
pt_file=['S151_m2_proximaltibia','S152_m3_proximaltibia', 'S152_m4_proximaltibia']


def merge_sph(fname,save):
	fw=open('plothist/'+save+'/spherical_coordinate.txt','w')
	fw1=open('plothist/'+save+'/spherical_coordinate_not_in_connected_component.txt','w')
	for i in range(len(fname)):
		f=open(fname[i]+'/'+'spherical_coordinate.txt','r')
		for line in f:
			l=line.split()
			l[5]=str(i+1)+l[5]
			l[6]=str(i+1)+l[6]
			#l=line[0:-1]+'\t'+str(i+1)+'\n'
			if int(l[1])==1:
				fw.write(l[0]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\n')
			else:
				fw1.write(l[0]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\n')	
	
def merge_deg(fname,save):
	fw=open('plothist/'+save+'/degree_of_clusters.txt','w')
	for i in range(len(fname)):
		f=open(fname[i]+'/'+'degree_of_clusters.txt','r')
		for line in f:
			fw.write(line)
	
		



merge_sph(df_file,'DF')
merge_sph(pt_file,'PT')


merge_deg(df_file,'DF')
merge_deg(pt_file,'PT')


