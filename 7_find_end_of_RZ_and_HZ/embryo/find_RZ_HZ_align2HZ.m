clear all 


allpath={
'data/Nuclei_and_CellsE185_S153_m7_distalfemur/',
'data/Nuclei_and_CellsE185_S153_m7_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m3_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m3_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m4_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m4_proximaltibia/',
}; 


Vdf=[];
for gi=[1,3,5]
    path=allpath{gi};
    s=strsplit(path,'Nuclei_and_CellsE185_');
     out2=strcat('MakeListNucleiLabelled/',s{2});    
     a=load([out2,'AllFeaturesSave_.mat']);
     Vdf=[Vdf;a.normalized_volume];
end

Vpt=[];
for gi=[2,4,6]
    path=allpath{gi};
    s=strsplit(path,'Nuclei_and_CellsE185_');
     out2=strcat('MakeListNucleiLabelled/',s{2});    
     a=load([out2,'AllFeaturesSave_.mat']);
     Vpt=[Vpt;a.normalized_volume];
end

Vdf=Vdf(find(Vdf<10));
Vpt=Vpt(find(Vpt<10));


subplot(1,2,1)
hist(Vdf,30);
subplot(1,2,2)
hist(Vpt,30)





for gi=5%:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        out1=strcat('MakeListNuclei_unlabelled/',s{2});    
        out2=strcat('MakeListNucleiLabelled/',s{2});    
        out3=strcat('MakeListClustersMask/',s{2});  
        a1=load([out1,'centroid_and_surface_nuclei.mat']);
        a2=load([out2,'centroid_and_surface_nuclei.mat']);
        a3=load([out3,'centroid_and_surface_nuclei.mat']);
        
        c1=a1.centroid - mean(a1.centroid);
        c2=a2.centroid - mean(a2.centroid);
        c3=a3.centroid - mean(a3.centroid);
        
        merge_centroid=[a1.centroid;a2.centroid;a3.centroid];
        merge_index1=(1:length(a1.centroid));  pp=length(a1.centroid);
        merge_index2=1+pp:pp+length(a2.centroid); pp=length(a1.centroid)+length(a2.centroid);
        merge_index3=1+pp:pp+length(a3.centroid); 
        
        merge_tileid=[a1.unique_tileid;a2.unique_tileid;a3.unique_tileid];
        section=unique(merge_tileid(:,3));
        shift_ori= merge_centroid-mean(merge_centroid);

        for se=1:length(section)
            id=find(merge_tileid(:,3)==section(se));   ankit(se,:)=[section(se),min(shift_ori(id,3)) max(shift_ori(id,3))];
        end
        
        maximum_shift=max(ankit(:,3));
        
        shift_coordinate=shift_ori;
        
        for se=1:length(section)
            id=find(merge_tileid(:,3)==section(se));   
            shift_coordinate(id,3)=shift_ori(id,3)+ maximum_shift - ankit(se,3);
            ankur(se,:)=[section(se),min(shift_coordinate(id,3)) max(shift_coordinate(id,3))];
        end
        
        ankur

        h=figure;
        plot3(shift_coordinate(merge_index1,1),shift_coordinate(merge_index1,2),shift_coordinate(merge_index1,3),'b.','markersize',3);
        hold on 
        plot3(shift_coordinate(merge_index2,1),shift_coordinate(merge_index2,2),shift_coordinate(merge_index2,3),'r.','markersize',3);
        plot3(shift_coordinate(merge_index3,1),shift_coordinate(merge_index3,2),shift_coordinate(merge_index3,3),'k.','markersize',3);
        xlabel('x');ylabel('y');zlabel('z')

        sname=s{2}(1:strlength(s{2})-1);
        saveas(h,['matlabFigures/', sname,'.fig'])
        saveas(h,['matlabFigures/',sname,'.png']);
        close all 
        

end
        