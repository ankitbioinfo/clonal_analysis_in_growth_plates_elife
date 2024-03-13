clear all 


dir='Nuclei_and_Cells_DT_S17_m2_wt/';

a=load([dir,'all_cells.mat']);
b=load([dir,'all_cells_nuclei.mat']);

cells=a.all_cells(:,5:7); %this is centroid 
nuc=b.all_cells_nuclei(:,5:7); %this is centroid 




file=load('cell_nuc_pairs.dat');
cell_pairid=file(:,1);
nuc_pairid=file(:,2);

%centroid=cells(cell_pairid,:); %use cell position as input 

[~,LCC,clcc,nlcc]=readClusterFile(['./RandomSampling/ulnaRealization/Random_cluster.dat'],cell_pairid,nuc_pairid); 


c=cells(clcc{5},:);
plot3(c(:,1),c(:,2),c(:,3),'bo');
hold on 
c=nuc(nlcc{5},:);
plot3(c(:,1),c(:,2),c(:,3),'r*');


mylcc=clcc;
cent=cells;

[length(clcc),length(nlcc)]

count=1;
for i=1:length(mylcc)
    if length(mylcc{i})>2
        clustermem=mylcc{i}';
        [CluAvgDeg(i,1), ColAvgDeg(i,1),largestDistance(i,1),stepsize,deg,edges,AZ_EL_R]=coordination_number(cent,clustermem);
        AllEdgesOFClusterSaved{count}=AZ_EL_R;
        count=count+1;
    end
end
 fid=fopen(['cell/spherical_coordinate.txt'],'w');
 for i=1:length(AllEdgesOFClusterSaved)
     E=AllEdgesOFClusterSaved{i};
     for j=1:size(E,1)
           fprintf(fid,'%d\t%d\t%0.3f\t%0.3f\t%0.3f\t%d\t%d\n',i,E(j,1),E(j,2),E(j,3),E(j,4),E(j,5),E(j,6));
     end
 end
 fclose(fid);
 clear AllEdgesOFClusterSaved;
 
 
 mylcc=nlcc;
cent=nuc;
 count=1;
for i=1:length(mylcc)
    if length(mylcc{i})>2
        clustermem=mylcc{i}';
        [CluAvgDeg(i,1), ColAvgDeg(i,1),largestDistance(i,1),stepsize,deg,edges,AZ_EL_R]=coordination_number(cent,clustermem);
        AllEdgesOFClusterSaved{count}=AZ_EL_R;
        count=count+1;
    end
end
 fid=fopen(['nuc/spherical_coordinate.txt'],'w');
 for i=1:length(AllEdgesOFClusterSaved)
     E=AllEdgesOFClusterSaved{i};
     for j=1:size(E,1)
           fprintf(fid,'%d\t%d\t%0.3f\t%0.3f\t%0.3f\t%d\t%d\n',i,E(j,1),E(j,2),E(j,3),E(j,4),E(j,5),E(j,6));
     end
 end
 fclose(fid);
 clear AllEdgesOFClusterSaved;
        
        

function [avgdeg, avgdegExpectedColumn,largestDistance,averageStepSize,deg,E,sph]=coordination_number(centroid,nodeid)
%function [avgdeg,diameter,averageStepSize]=coordination_number(centroids,edges)

%              id=unique(edges(:));
%              for j=1:length(id)
%                  new2old(j)=id(j);
%                  old2new(id(j))=j;
%              end
%                bondDist=[];
%                for i=1:size(edges,1)
%                  newedges(i,1)=old2new(edges(i,1));
%                  newedges(i,2)=old2new(edges(i,2));
%                  bondDist=[bondDist;pdist(centroids(edges(i,:),:))];
%                end
            cent=centroid(nodeid,:);
            n=size(cent,1);
            dist=[];
            myedge=[];
            for i=1:n
                for j=i+1:n
                    dist=[dist,pdist(cent([i,j],:))];
                    myedge=[myedge;[i,j]];
                end
            end
            
            [sa,sb]=sort(dist);
            %[size(myedge), size(dist)]
            
            start=1;
            for i=2:length(sa)
                nodes=unique(myedge(sb(1:i),:));
                % [i,n,length(nodes)]
                if n==length(nodes)
                    start=i;
                    break 
                end
            end
            
          
            
            flag=true;
            while flag
                E=myedge(sb(1:start),:);
                G=graph(E(:,1),E(:,2));
                bins=conncomp(G);
                % number of connected components 
                nocomp=unique(bins);
                clear numberOfObjectsInConnectedComponents
                for i=1:length(nocomp)
                    numberOfObjectsInConnectedComponents(i)=sum(nocomp(i)==bins);
                end

                %[start,numberOfObjectsInConnectedComponents]
                if length(nocomp)==1
                    flag=false;
                    TotalE=start; 
                end
                start=start+1;
            end
               
            
            E=myedge(sb(1:TotalE),:);
            averageStepSize=mean(sa(1:TotalE));

            largestDistance=sa(TotalE);
            G=graph(E(:,1),E(:,2));
            deg=G.degree();
            avgdegofcompletegraph=length(deg)-1;
            avgdeg=mean(deg);
            avgdegExpectedColumn=2*(n-1)/n;
            %avgdeg= mean(deg)/degreeExpectedColumn; 
            %[length(nodes),length(deg)]
%             
%             d=distances(G);
%             diameter=max(d(:));

%             c=cent-mean(cent);
%             [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
%             sph=[az,el,r];
            
            count=1;
            for i=1:size(E,1)
                c=cent(E(i,:),:);
                c=c-mean(c);
                [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
                flag=[1;1];
                sph(count:count+1,:)=[flag,az,el,r,nodeid(E(i,[1,2]))',nodeid(E(i,[2,1]))'];
                count=count+2;
            end    
                
             for i=1:size(myedge,1)
                c=cent(myedge(i,:),:);
                c=c-mean(c);
                [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
                flag=[0;0];
                sph(count:count+1,:)=[flag,az,el,r,nodeid(myedge(i,[1,2]))',nodeid(myedge(i,[2,1]))'];
                count=count+2;
            end    
            

end



 function [LCC,LCC1,cell_lcc,nuc_lcc]=readClusterFile(name, cellid, nucid) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
         while ischar(tline)
            line= split(tline);
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
        
        for i=1:length(LCC1)
            myid=LCC1{i};
            cell_lcc{i}=cellid(myid);
            nuc_lcc{i}=nucid(myid);
        end
        
 end

