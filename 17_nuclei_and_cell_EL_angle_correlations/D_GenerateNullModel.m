
clear all 
close all

warning('off','all')

randomSampling='RandomSampling/';
if ~exist([randomSampling],'dir')
     mkdir([randomSampling]);
end



allpath={
'data/Nuclei_and_CellsE185_S153_m7_distalfemur/',
'data/Nuclei_and_CellsE185_S153_m7_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m3_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m3_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m4_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m4_proximaltibia/',
}; 

RZ_HZ_height={ [-429,632]};



for gi=1%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        s={'ulna','ulna'};
        outputpath_unlabelled=strcat('MakeListNuclei_unlabelled/',s{2});
        outputpath_labelled=strcat('MakeListNucleiLabelled/',s{2});
        %outputpath_mask=strcat('MakeListClustersMask/',s{2});
        rand_dir_output=strcat(randomSampling,s{2});
        
        if ~exist([rand_dir_output],'dir')
            mkdir([rand_dir_output]);
        end
              
        dir2=strcat(rand_dir_output,'Realization','/');
        if ~exist([dir2],'dir')
            mkdir([dir2]);
        end
end        


%load([outputpath_unlabelled,'centroid_and_surface_nuclei.mat'],'centroid');
file=load('cell_nuc_pairs.dat');
dir='Nuclei_and_Cells_DT_S17_m2_wt/';
a=load([dir,'all_cells.mat']);
cells=a.all_cells(:,5:7);
centroid=cells(file(:,1),:); %use cell position as input 
%cell_cutoff= load('../Tibia/MakeListColumnarStructurePrediction/cells_cutoff.dat');



radius=load(['ParameterForRandomModel_basedon_nuc_.dat']);
%randomCellCutoff=cell_cutoff(:,2);
[~,LCC]=readClusterFile(['degree_of_clusters.txt']); 
disp('LCC')
length(LCC)

% Part 1 
c=centroid-mean(centroid);
%bonetype=1;
PDaxis=3; DV=1; ML=2;
% if (bonetype==3)|(bonetype==1)
%        c=[-c(:,1),c(:,2:3)];
% else
%        c=[c(:,1:2),c(:,3)];
% end

plot3(c(:,1),c(:,2),c(:,3),'b.');
xlabel('X'); ylabel('Y'); zlabel('Z'); 
RZ=min(c(:,PDaxis));
HZ=max(c(:,PDaxis));
axis image


RZ=RZ_HZ_height{gi}(1);HZ=RZ_HZ_height{gi}(2);
[RZ, HZ]


% Part 2 
gaussian=fittype(@(mu,sigma,x)  ( 1/(sigma*sqrt(2*pi))  *exp( -0.5*(  ((x-mu)./sigma).^2 ) )    ));

[y1,x1]=histnorm(c(:,1),15);
[y2,x2]=histnorm(c(:,2),15);
[y3,x3]=histnorm(c(:,3),15);

X=[x1',x2',x3'];
Y=[y1',y2',y3'];

tname={'X','Y','Z'};
for i=1:3
    xr=X(:,i); yr=Y(:,i);
    [f,stat]=fit(xr,Y(:,i),gaussian,'startpoint',[mean(c(:,i)),std(c(:,i))]);
    fit_mu_std{i}=f;
end

% h=figure;
% mycolor={'r.-','g.-','b.-'};
% for i=1:3
%     p(i)=plot(X(:,i),Y(:,i),mycolor{i});hold on 
%     xr=X(:,i); f= fit_mu_std{i};
%     plot(xr,f(xr),'k-')
%     leg{i}=strcat(tname{i},':',', fit(\mu=',sprintf('%0.3f',f.mu),', \sigma=',sprintf('%0.5f',f.sigma),')');
% end
% 
% legend(p,leg,'location','north')
% xlabel('space')
% ylabel('probability')

clear p 
clear leg




% Part 3 
mu_std{1}=[fit_mu_std{1}.mu, fit_mu_std{1}.sigma];
mu_std{2}=[fit_mu_std{2}.mu, fit_mu_std{2}.sigma];
mu_std{3}=[fit_mu_std{3}.mu, fit_mu_std{3}.sigma];
 


 
[yr,xr]=histnorm(radius(:,1),10);
[f,stat]=fit(xr',yr',gaussian,'startpoint',[mean(radius(:,1)),std(radius(:,1))]);
% h=figure;
% p(2)=plot(xr,f(xr),'k-');leg{1}='data';
% hold on 
% p(1)=plot(xr,yr,'r.-'); 
% leg{2}=strcat('fit(\mu=',sprintf('%0.3f',f.mu),', \sigma=',sprintf('%0.5f',f.sigma),')');
% legend(p,leg,'location','north')




RR=[];
for iteration=1:2
    fid=fopen([dir2,'Random_cluster_',num2str(iteration),'.dat'],'w');
       
    dir1=strcat(dir2,'Iteration_',num2str(iteration),'/');
    if ~exist([dir1],'dir')
       mkdir([dir1]);
    end   
    %disp('length of radius')
    %length(radius)
    for k=1:length(radius)
        [iteration,k]
        flag=1;
        while flag
            %R(k,1)= f.mu + f.sigma.*randn(1);
            %R(k,1)=radius(k,1)+2*0.05*rand-0.05;  % according to 21 bin length 
            %R(k,1)=radius(k,1)+2*0.1429*rand-0.1429;  % according to 8 bin length
            R(k,1)=radius(k,1)+2*0.25*rand-0.25;  % according to 5 bin length

            
            if (R(k,1)>=0)&(R(k,1)<=1)
                 Z=RZ+(HZ-RZ)*R(k,1);     
                 X=mu_std{DV}(1) + mu_std{DV}(2)*randn(1);
                 Y=mu_std{ML}(1) + mu_std{ML}(2)*randn(1);
                 size_of_cluster=length(LCC{k});
                 random_radius_factor= 0.1*size_of_cluster*rand*radius(i,2); 
                 cellid=Search_cell_in_approximity(X,Y,Z,random_radius_factor,c,DV,ML,PDaxis);
                 %[k, X, Y, Z, length(cellid), size_of_cluster]
                  if (length(cellid)-size_of_cluster)>=0
                        
                        binindex=which_bin(R(k,1));
                        cutoff=randi(20);
                        %cutoff=randomCellCutoff(binindex);
                      
                        [cid,failure1,edges]=find_small_random_cluster(cellid,c,size_of_cluster,cutoff);
                        if failure1==0
                           %[iteration,k,length(cellid),length(cid),size_of_cluster,failure]
                             [cid,nedges,failure2]=select_few_connected_cells_from_pool(edges,cid,size_of_cluster);
                             if failure2==0
                                    flag=0;
                                    for j=1:length(cid)
                                        fprintf(fid,'%d\t',cid(j));
                                    end
                                    fprintf(fid,'\n'); 
                                    %[size(edges,1),size(nedges,1)]
                                    clusterEdges{k}=nedges;
                                    nocomp=number_of_connected_componenets(nedges);
                                    if nocomp~=1                                    
                                           
                                                k
                                                cid
                                                edges
                                                size_of_cluster
                                                pause
                                                                           
                                    end
                                    
                             end
                                                        
                        end
                         
                  end
            end
        end
    end
        
    fclose(fid);    
    
    length(clusterEdges)
    
    
    for k=1:length(clusterEdges)
        nedges=clusterEdges{k};
        fid1=fopen(strcat(dir1,'Edges_',num2str(k),'.dat'),'w');
        for i=1:size(nedges,1)
            fprintf(fid1,'%d\t%d\n',nedges(i,1),nedges(i,2));
        end
        fclose(fid1);
    end
    
    
    
    if iteration<10
        RR=[RR;R];
    end
   
end


%save('RanomClusterEdges.mat','clusterEdges')


%  [yr,xr]=histnorm(RR,20);
%  plot(xr,yr,'bo-');


%binsize=51;zone=linspace(RZ,HZ,binsize);Xbin=linspace(0,1,binsize-1);


function binindex=which_bin(data)
    binindex=-1;
    interval=linspace(0,1,11);
    for i=1:length(interval)-1
        if ((data>=interval(i))&(data<interval(i+1)))
                binindex=i;
        end
    end
end


function [cid,nedges,failure]=select_few_connected_cells_from_pool(edges,given_cid,size_of_cluster)
     failure=0;
     [~,ia]=unique(edges,'rows');
     edges=edges(ia,:);
     v=1:length(given_cid);
     cid=given_cid(v(1:size_of_cluster));nedges=search_edges(edges,cid);
     cid=unique(nedges(:));
     
     if length(nedges)==0
         failure=1; nocomp=0;
     else
         nocomp=number_of_connected_componenets(nedges);
     end
     
     while (nocomp~=1)
           nv=v(randperm(length(v)));
           cid=given_cid( nv(1:size_of_cluster));nedges=search_edges(edges,cid);
           cid=unique(nedges(:));
           if length(nedges)==0
                failure=1;
           else
                nocomp=number_of_connected_componenets(nedges);
           end
     end
     
     if length(cid)~=size_of_cluster
         failure=1;
     end
     
end
     
                          
                        
function myedges=search_edges(edges,vertices)
%     edges
%     vertices 
    [~,ia]=unique(edges,'rows');
    edges=edges(ia,:);
	count=1;
    myedges=[];
	for j=1: size(edges,1)
		flag1=0;
		flag2=0;
		for i=1:length(vertices)
			if (edges(j,1)==vertices(i))
				flag1=1;
            end
			if (edges(j,2)==vertices(i))
				flag2=1;
            end
            if (flag1+flag2)==2
                myedges(count,:)=edges(j,:);
                count=count+1;
            end
        end
    end
    
    [~,ia]=unique(myedges,'rows');
    myedges=myedges(ia,:);
    
 end




function cellid=Search_cell_in_approximity(X,Y,Z,radius,c,DV,ML,PDaxis)
  distance=(c(:,DV) -  X).^2  +  (Y - c(:,ML)).^2 + (Z - c(:,PDaxis)).^2;
  cellid=find(distance< (radius^2));
end

function [cid,failure,edges]=find_small_random_cluster(cellid,centroid,size_of_cluster,cutoff)

          failure=0;
          edges=[];
          cid=[];
          maximum_iteration=500;
          c=centroid(cellid,:);
          count=1;
          dist=[];
          for i=1:length(cellid)
              for j=i+1:length(cellid)
                  d=pdist(c([i,j],:));
                  if d<cutoff
                      dist(count,:)=[d,cellid(i),cellid(j)];
                      count=count+1;
                  end
              end
          end
          
          if length(dist)>0
                  tempcell=dist(:,2:3);
                  if length(unique(tempcell(:)))<size_of_cluster
                       failure=1;
                       
                  else     

                          [sa,sb]=sort(dist(:,1));

                          sorted_dist=dist(sb,2:3);
                          cid=unique(dist(sb(1),[2:3]));
                          myindex=find_matching_edges(sorted_dist,cid);
                          edges=sorted_dist(myindex,:);


                          if size_of_cluster~=length(cid)
                              flag=1;
                              count=0;
                              while flag==1
                                      count=count+1;
                                      myindex=find_matching_edges(sorted_dist,cid);
                                      %[size(sorted_dist),length(cid),length(myindex)]
                                      edges=sorted_dist(myindex,:);
                                      nocomp=number_of_connected_componenets(edges);
                                      cid=unique(edges(:));
                                      if (length(cid)>=size_of_cluster)
                                          %[size_of_cluster,length(cid),length(cellid),nocomp]
                                          flag=0;
                                          %pause
                                      end

                                      if count>maximum_iteration
                                          failure=1;
                                          flag=0;
                                      end

                              end
                          end 
                  end
          else
                 failure=1;
          end
          
end


function finalIndex=find_matching_edges(edges,cid)
          finalIndex=[];
          for i=1:length(cid)
              myindex= find((cid(i)==edges(:,1))|(cid(i)==edges(:,2)));
              finalIndex=[finalIndex;myindex];
          end
          finalIndex=unique(finalIndex);
end
              
                  
             

function LCC=LargestConnectedComponents(edges)

        [~,ia]=unique(edges,'rows');
        edges=edges(ia,:);
        
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



function nocomp=number_of_connected_componenets(edges)
         [~,ia]=unique(edges,'rows');
         edges=edges(ia,:);

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
        nocomp=length(unique(bins));
end

          
          
          

 function [LCC,LCC1]=readClusterFile(name) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
         while ischar(tline)
            line= split(tline,',');
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




