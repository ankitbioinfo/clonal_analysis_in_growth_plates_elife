


df_file=['S153_m7_distalfemur','S154_m3_distalfemur','S154_m4_distalfemur']
pt_file=['S153_m7_proximaltibia','S154_m3_proximaltibia', 'S154_m4_proximaltibia']

fname=df_file
fw=open('ParameterForRandomModel_basedon_nuc_.dat','w')
fw2=open('degree_of_clusters.txt','w')
for i in range(len(fname)):
	f=open('./MakeListNucleiLabelled/'+fname[i]+'/ParameterForRandomModel_basedon_nuc_.dat')
	for line in f:
		fw.write(line)
	f=open('./MakeListNucleiLabelled/'+fname[i]+'/degree_of_clusters.txt')
	for line in f:
		fw2.write(line)	
