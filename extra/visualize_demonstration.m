clear all 


allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 


path=allpath{6};
s=strsplit(path,'Nuclei_and_CellsP40_');
outputpath=strcat('../proximalTibia_DistalFemur/MakeListNucleiLabelled/',s{2});

a=load([outputpath,'centroid_and_surface_nuclei.mat']);


id=[7962
       8027
       8033
       8066
       8078
       8105
       8150
       8157
       8172
       8184
       8194
       8274
       8280
       8281
       8286
       8317
       8324];
   
 id=   [3118,
       3137,
       4866];  %EL > 80 and AZ >0 and AZ<10 
   
  id=[  3815
       3925
       4026]; %EL > 80 and AZ >0 and AZ<10 
   
  id=[    18428
      18451
      18452
      18463]; %EL > 80 and AZ >80 and AZ<90 
  
  id=[    991
       1004
       1016];  %EL > 80 and AZ >80 and AZ<90 
   
  id=[
       6763
       6796
       6828]; %EL > 80 and AZ >0 and AZ<10 
   
  id=[
       8194
       8286
       8324]; 
   
  h=figure; 
  mycolor={'c.','r.','y.'}; 
  for i=1:length(id)
      c=a.nuc{id(i)};
      ct=mycolor{a.unique_tileid(id(i),1)};
      plot3(c(:,1),c(:,2),c(:,3),ct,'markersize',4);   
      hold on 
      cent=mean(c);
      text(cent(1),cent(2),cent(3),num2str(id(i)),'fontsize',8);
  end
  
  
xlabel('x')
ylabel('y')
zlabel('z') 
saveas(h,['cluster_AZ_EL_','.png']); 
saveas(h,['cluster_AZ_EL_','.fig']); 

