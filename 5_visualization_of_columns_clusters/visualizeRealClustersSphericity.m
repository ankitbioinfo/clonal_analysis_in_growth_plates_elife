

clear all 

clear all 
warning('off','all')


allpath={
'data/Nuclei_and_CellsE185_S153_m7_distalfemur/',
'data/Nuclei_and_CellsE185_S153_m7_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m3_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m3_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m4_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m4_proximaltibia/',
}; 


for gi=1%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        outputpath_labelled=strcat('MakeListNucleiLabelled/',s{2});
        colmask=load(['MakeListClustersMask/',s{2},'centroid_and_surface_nuclei.mat']);
        %outputpath_mask=strcat('MakeListClustersMask/',s{2});
        %outputpath_unlabelled=strcat('MakeListNuclei_unlabelled/',s{2});
        %rand_dir_output=strcat('RandomSampling/',s{2});

end        

%radius=load([outputpath_labelled,'ParameterForRandomModel_basedon_nuc_.dat']);
nuclei=load([outputpath_labelled,'centroid_and_surface_nuclei.mat']);


a=load([outputpath_labelled,'AllFeaturesSave_.mat']);
sphericity=a.sphericity; 


savepath=strcat(outputpath_labelled,'Visualization','/');
if ~exist([savepath],'dir')
    mkdir([savepath]);
end

[LCC,allid]=readClusterFile(strcat(outputpath_labelled,'Cluster_NoDuplication.dat'));

mycolor={'c.','r.','y.'};
for i=1:100%length(sphericity)
    id=LCC{i};
    h=figure;
    for j=1:length(id)
         v=nuclei.nuc{id(j)};
         ct=mycolor{nuclei.unique_tileid(id(j),1)};
         plot3(v(:,1),v(:,2),v(:,3),ct);
         hold on 
    end
    
    sp=round(100*sphericity(i));
    
    title(['# of cells ',num2str(length(id)) , ' , sphericity ', sprintf('%0.2f',sphericity(i)) ]);
    saveas(h,[savepath,'cluster',num2str(sp),'_',num2str(i),'.png']);
    saveas(h,[savepath,'cluster',num2str(sp),'_',num2str(i),'.fig']);
    close all 
end


 function [LCC,allCellId]=readClusterFile(name) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
        allCellId=[];
         while ischar(tline)
            line= split(tline);
             if length(line)>2
                for j=1:length(line)-1
                    LCC{count}(j)=str2num(line{j});
                    allCellId=[allCellId;LCC{count}(j)];
                end
             end    
%             
%              if length(line)>2
%                 for j=1:length(line)-1
%                     LCC1{count}(j)=str2num(line{j});
%                 end
%              end    
             
            
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
 end
    
    