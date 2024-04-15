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

RZ_HZ_height={ [-259,486], [-316,613], [-339,511],[-277,516],[-293,654],[-259,500]};

outputVisualize='merge_clones/';
if ~exist([outputVisualize],'dir')
     mkdir([outputVisualize]);
end

for gi=6%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        inputpath=strcat('../MakeListNucleiLabelled/',s{2});
        inputpath_mask=strcat('../MakeListClustersMask/',s{2});
        outputpath=strcat(outputVisualize,s{2});
end

if ~exist([outputpath],'dir')
     mkdir([outputpath]);
end


load([inputpath,'centroid_and_surface_nuclei.mat']);
sarah_mask=load([inputpath_mask,'centroid_and_surface_nuclei.mat']);
sarah_mask_id=load([inputpath,'ClusterMask_NoDuplication.dat']);

[~,LCC]=readClusterFile([inputpath,'Cluster_NoDuplication.dat']); 

savepath='figures/';

edges=[];
for j=1:length(LCC)
    id1=LCC{j};
    
    c1=centroid(id1,:);
    
    %[cell1,comb1]=FindShapeCluster(id1,nuc);
    
    
    
     for k=1:length(LCC)
        id2=LCC{k};
        c2=centroid(id2,:);
        d=pdist2(c1,c2);
        neighborCluster(k)=min(d(:));
        if k==j
            neighborCluster(k)=100000000;
        end
     
        
%        [cell2,comb2]=FindShapeCluster(id2,nuc);
%        combined=[cell1;cell2];   
%        [~,combV]=convhull(combined);
%         volumeFraction= combV/(comb1+comb2);
%         if volumeFraction<1
%             h=figure;
%             
%             plot3(cell1(:,1),cell1(:,2),cell1(:,3),'b.');
%             hold on 
%             plot3(cell2(:,1),cell2(:,2),cell2(:,3),'r.');
%              saveas(h,[savepath,'cluster',num2str(j),'_',num2str(k),'.png']);
%              %saveas(h,[savepath,'cluster',num2str(myid),'.fig']);
%             close all 
%         end

    
     end
     
      [sa,sb]=min(neighborCluster);
      if sa<20
          edges=[edges;[j,sb]];
      end
     ankit(j,:)=[sa,sb];
end

mergeLCC=LargestConnectedComponents(edges);
a=unique(edges(:));
b=1:length(LCC);
c=setdiff(b,a);

count=1+length(mergeLCC);
for i=1:length(c)
    mergeLCC{count}=c(i);
    count=count+1;
end

[length(LCC), length(mergeLCC)]


fid=fopen([ outputpath,'Cluster.dat'],'w');
for i=1:length(mergeLCC)
    nodes=[];
    for j=1:length(mergeLCC{i})
        id=mergeLCC{i}(j);
        nodes=[nodes,LCC{id}];
    end
    
    id_notgood=nodes;
    
    Repeat_centroids=centroid(id_notgood,:);
    [~,firstgood]=unique(Repeat_centroids,'rows');
     id=id_notgood(firstgood);
     Repeat_centroids=centroid(id,:);
     Repeat_volume=celvolume(id,1);
     clear Repeat_surfaces;
     for j=1:length(id)
         Repeat_surfaces{j}=nuc{id(j)};
     end
            
     ia=RemoveBadCells(Repeat_centroids,Repeat_volume,Repeat_surfaces);
     LCC_nuc{i}=id(ia); 
    
end
    

LCC_nuc= sorted_cluster_from_larger_to_smaller(LCC_nuc);

for i=1:length(LCC_nuc)
    for j=1:length(LCC_nuc{i})
        fprintf(fid,'%d ',LCC_nuc{i}(j));
    end
    fprintf(fid,'\n');
end
fclose(fid);




function [LCC]=sorted_cluster_from_larger_to_smaller(C)
    LCC={};
    count=1;
    for i=1:length(C)
         clulen(count,1)=length(C{i});
         count=count+1;
    end
    [sa,sb]=sort(clulen,'descend');
    
    for i=1:length(sb)
        LCC{i} = C{ sb(i)};
    end
    %newmask=mask(sb);
end


function [ClusterShape,V]=FindShapeCluster(id,nuc)
     mergeVol=0;
     ClusterShape=[];
     for i=1:length(id)
            ClusterShape=[  ClusterShape;nuc{id(i)}];
            %mergeVol=mergeVol+celvolume(id(i));
     end
     %combinedCellVolumeInCluster(j,1)=mergeVol;
     [K,V]=convhull(ClusterShape);
end



function goodcellindex=RemoveBadCells(centroid, volume,surfaces)
        if size(centroid,1)>4
            [~,firstgood]=unique(centroid,'rows');
                %tic
                %[length(centroid),length(firstgood)]
            % spacing=[1, 1, 1];
            % delta_spacing=2;
            % if the user didn't define 'exclude_boundary_points' we set it to false:
            if ~exist('exclude_boundary_points', 'var')
                exclude_boundary_points = false;
            end
            % calculating the triangulation and the volume of each triangle:
            N=centroid(firstgood,:);

            TRI = delaunay(N(:,1), N(:,2), N(:,3));
            clear neighbor
            edges=[];
            count=1;
            for i = 1 : size(N,1)
                temp=[];
                for j=1:size(TRI,1)
                    for k=1:size(TRI,2)
                        if TRI(j,k)==i
                            temp=[temp,TRI(j,:)];
                        end
                    end
                end


                neighborList=setdiff(unique(temp),i);
                %neighbor(i,1)=length(neighborList{i});
                ids= neighborList;

                for j=1:length(ids)
                    a=min(ids(j),i);
                    b=max(ids(j),i);
                    edges(count,:)=[firstgood(a),firstgood(b)];
                    count=count+1;
                end

            end

            [~,ia]=unique(edges,'rows');
            edges=edges(ia,:);
        else
            edges=[];
            n=size(centroid,1);
            firstgood=1:n;
            for i=1:n 
                for j=i+1:n
                    edges=[edges;[i,j]];
                end
            end
        end    
            
%             edges
%              length(centroid)

             badcell=[];
             for j=1:size(edges,1)
                 dist=pdist(centroid(edges(j,:),:));
                 if dist<10
                         cell1=surfaces{edges(j,1)}; actualVol(1)=volume( edges(j,1)   );      
                         cell2=surfaces{edges(j,2)}; actualVol(2)=volume( edges(j,2)   );  
                         combined=[cell1;cell2];
                         [~,combV]=convhull(combined);
                         [~,comb1]=convhull(cell1);
                         [~,comb2]=convhull(cell2);

                         volumeFraction= combV/(comb1+comb2);
                         if volumeFraction<1
                                    [sa,sb]=min(actualVol);
                                    if sb==1
                                        badid=edges(j,1);
                                    end
                                    if sb==2
                                        badid=edges(j,2);
                                    end
                                    badcell=[badcell;badid];
                         end
                 end
             end    

              goodcellindex=setdiff(firstgood,badcell);
              %toc
              %length(goodcellindex)  
end





 function [LCC,LCC1]=readClusterFile(name) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
         while ischar(tline)
            line= split(tline,' ');
             if length(line)>3
                for j=1:length(line)-1
                    LCC{count}(j)=str2num(line{j});
                end
             end    
            
             if length(line)>2
                for j=1:length(line)-1
                    LCC1{count}(j)=str2num(line{j});
                end
             end    
             
            
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
 end
 
 
 

function LCC=LargestConnectedComponents(edges)
        cellIds=unique(edges(:));
        for j=1:length(cellIds)
            old2new(cellIds(j),1)=j;
            new2old(j,1)=cellIds(j);    
        end
        [length(old2new),length(new2old),length(cellIds)];

        for i=1:size(edges,1)
            for j=1:2 
                newedgename(i,j)= old2new(edges(i,j));
            end
        end
        G=graph(newedgename(:,1),newedgename(:,2));
        bins=conncomp(G);
        % number of connected components 
        nocomp=unique(bins);
        %disp(['# of connected components  ', num2str(length(nocomp))]);
        for i=1:length(nocomp)
            numberOfObjectsInConnectedComponents(i)=sum(nocomp(i)==bins);
        end
        
        
        [sa,sb]=sort(numberOfObjectsInConnectedComponents,'descend');
        
        index=1;
        for i=1:length(sa)
            if sa(i)>1
                LCCIds=find(bins==nocomp(sb(i)));
                LCC{index}=new2old(LCCIds);
                index=index+1;
            end
        end
end
