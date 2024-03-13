

clear all
allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 


outputVisualize='visualize_Eli/';
if ~exist([outputVisualize],'dir')
     mkdir([outputVisualize]);
end




for gi=6   %2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        outputpath=strcat('MakeListNucleiLabelled/',s{2});
        
        colmask=load(['MakeListClustersMask/',s{2},'centroid_and_surface_nuclei.mat']);
        %tileid=unique(colmask.unique_tileid(:,3));
        a=load([outputpath,'centroid_and_surface_nuclei.mat']);
        %cloneid=unique(colmask.unique_tileid(:,1));
        
end        

savepath=strcat(outputVisualize,s{2});    
        if ~exist([savepath],'dir')
             mkdir([savepath]);
end



 
%[LCC,LCC1]=readClusterFile(['ClusterGreaterThan100Micron1.dat']);
[LCC,LCC1]=readClusterFile([outputpath,'Cluster_good.dat']);
colmaskid=load([outputpath,'ClusterMask_good.dat']);

%col_like_cells=[1496,1526,1561,1584,1614,1657,1696];


for i=1:length(LCC)
    id=LCC{i};
    cent=a.centroid(id,:);
    t=unique(a.unique_tileid(id,1));
    if length(t)>1
        for k=1:length(t)
            temp=find(t(k)==a.unique_tileid(id,1));
            ankur{i}{k}=id(temp);
        end
    end
    %clen(i,:)=max(cent)-min(cent);
end


for myid= 1:length(LCC)
id=LCC{myid};
posid=unique(a.unique_tileid(id,3));

h=figure;

mycolor={'c.','r.','y.'};
%yellow replaced by black 
ClusterShape=[];
mergeVol=0;

[E1,sph,deg,mycolor]=coordination_number(a.centroid,id);
connectEdg= id(E1); 

for j=1:size(connectEdg,1)
    E=connectEdg(j,:);
    plot3(a.centroid(E,1), a.centroid(E,2), a.centroid(E,3), 'K-','linewidth',2);
    hold on 
end
    
    
BoneGrowthAxis=[0; 0 ; 1];

for i=1:length(id)
    c=a.nuc{id(i)};
    %a.unique_tileid(id(i),:)
    
    ct=mycolor{i};
    %ct=mycolor{a.unique_tileid(id(i),1)};
    plot3(c(:,1),c(:,2),c(:,3),ct,'markersize',5);
   
    hold on 
    cent=a.centroid(id(i),:);
    %ankit(i,:)=[id(i),cent];
    %ankit(i,:)=max(cent)-min(cent);
   % text(cent(1),cent(2),cent(3),num2str(id(i)),'fontsize',8);
    
%     hold on 
     plot3(cent(1),cent(2),cent(3),'k*')
    
%     Repeat_centroids(i,:)=a.centroid(id(i),:);
%     Repeat_volume(i,:)=a.celvolume(id(i),:);
%     Repeat_surfaces{i}=c;
    
  %   ClusterShape=[  ClusterShape;c];
     %mergeVol=mergeVol+celvolume(LCC{j}(i));
    
end 



%height=max(ankit(:,4))-min(ankit(:,4))
%nocell=length(ankit);
%[max(ankit(:,4)),min(ankit(:,4))]

xlabel('x')
ylabel('y')
zlabel('z')
%title(['# of cells ',num2str(length(id)) , ', section id ',num2str(posid) , ' , cluster height ', num2str(clen(myid,3)) ])
%view(6,2); %67
%view(-10,-3); %357
%view(-96,0) %136
%[caz,cel]=view

%view(29.4,1.8) %clone 3 cluster 1 
%view(14.966,3.02)%clone 3 cluster 2 

c=colmask.nuc{colmaskid(myid)};
[clusterPCvector,~,latent]=pca(c);
AngleBetweenClusterPC1AndBone_PD=90-oangle(clusterPCvector(:,1),BoneGrowthAxis);
%plot3(c(:,1),c(:,2),c(:,3),'k.','markersize',0.75);

    [ center, radii, evecs,v, ~ ] = ellipsoid_fit_new( c );
     %draw fit
    mind = 0.5*min( c );
    maxd = 1.2*max( c );
    %nsteps = 200;     step = ( maxd - mind ) / nsteps;
    step=[1.9,1.9,1.9];
    n1=length(mind(1):step(1):maxd(1));
    n2=length(mind(2):step(2):maxd(2));
    n3=length(mind(3):step(3):maxd(3));
    nsteps=[n1,n2,n3];
    a1=linspace( mind(1) - step(1), maxd(1) + step(1), nsteps(1) );
    a2=linspace( mind(2) - step(2), maxd(2) + step(2), nsteps(2) );
    a3=linspace( mind(3) - step(3), maxd(3) + step(3), nsteps(3) ); 
    
    [ x, y, z ] = meshgrid(a1 ,a2 ,a3 );
    Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
              2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
              2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
    
  
    p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );
    set( p, 'FaceColor', 'm', 'EdgeColor', 'none' ); alpha(0.1);
    camlight;
     lighting phong;      






if AngleBetweenClusterPC1AndBone_PD>60
    savename=strcat('column_', num2str(ceil(AngleBetweenClusterPC1AndBone_PD)) );
    savename2= strcat(', column =', sprintf('%0.1f',AngleBetweenClusterPC1AndBone_PD));
else
    savename=strcat(' cluster_', num2str(ceil(AngleBetweenClusterPC1AndBone_PD)));
    savename2= strcat(', cluster =', sprintf('%0.1f',AngleBetweenClusterPC1AndBone_PD));

end

title(['# of cells ',num2str(length(id)) , ', section id ',num2str(posid) ,savename2 ])


view(-94.64,11.29)

saveas(h,[savepath,savename,'_id_',num2str(myid),'.png']);
saveas(h,[savepath,savename,'_id_',num2str(myid),'.fig']);
close all 

end




 function  angle= oangle(u,v)
      angle=atan2(norm(cross(u,v)),dot(u,v));
      if angle>(pi/2)
          angle=pi-angle;
      end
      factor=180/pi;
      angle=angle*factor;
 end

 %ia=RemoveBadCells(Repeat_centroids,Repeat_volume,Repeat_surfaces);


function goodcellindex=RemoveBadCells(centroid, volume,surfaces)
    tic
    [~,firstgood]=unique(centroid,'rows');
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

%size(edges)
[~,ia]=unique(edges,'rows');
edges=edges(ia,:);
%size(edges)

 badcell=[];
 for j=1:size(edges,1)
     edges(j,:)
     dist=pdist(centroid(edges(j,:),:));
     if dist<10
             cell1=surfaces{edges(j,1)}; actualVol(1)=volume( edges(j,1)   );      
             cell2=surfaces{edges(j,2)}; actualVol(2)=volume( edges(j,2)   );  
             combined=[cell1;cell2];
             [~,combV]=convhull(combined);
             [~,comb1]=convhull(cell1);
             [~,comb2]=convhull(cell2)

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
  toc
  %length(goodcellindex)  
  
end


function [E,sph,deg,mycolor]=coordination_number(centroid,nodeid)
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
            G=graph(E(:,1),E(:,2));
            deg=G.degree();
 
            count=1;
            for i=1:size(E,1)
                c=cent(E(i,:),:);
                c=c-mean(c);
                [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
                flag=[1;1];
                sph(count:count+1,:)=[flag,az,el*180/pi,r,(E(i,[1,2]))',(E(i,[2,1]))'];
                count=count+2;
            end    
            
            nuc_EL_Color = zeros(max(E(:)),1) ;
            degCount = zeros(max(E(:)),1) ;

            
            
             for i=1:size(sph,1)
                 if sph(i,3)>0
                     nuc_EL_Color(sph(i,5),1) =   nuc_EL_Color(sph(i,5),1) + sph(i,3);
                     nuc_EL_Color(sph(i,6),1) =   nuc_EL_Color(sph(i,6),1) + sph(i,3);
                     
                     degCount(sph(i,5),1) =   degCount(sph(i,5),1) + 1;
                     degCount(sph(i,6),1) =   degCount(sph(i,6),1) + 1;
                 end
             end
             
             
             avgEL = nuc_EL_Color./degCount;
             
             bininterval=[0,40,50,60,70,80,92];
             %bininterval=[0,10,20,30,40,50,92];
             %mycolorAll={'b','r','m','g','c','y'};
             mycolorAll={'y.','c.','g.','m.','r.','b.'};
             
             mycolor=cell(length(avgEL),1);
             
             
             for i=1:length(avgEL)
                 for j=1:length(bininterval)-1
                     if ((bininterval(j) <=  avgEL(i)) & (avgEL(i)<  bininterval(j+1)))
                         mycolor{i}=mycolorAll{j};
                     end
                 end
             end

                    
                
%             for i=1:size(myedge,1)
%                 c=cent(myedge(i,:),:);
%                 c=c-mean(c);
%                 [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
%                 flag=[0;0];
%                 sph(count:count+1,:)=[flag,az,el,r,nodeid(myedge(i,[1,2]))',nodeid(myedge(i,[2,1]))'];
%                 count=count+2;
%             end    
%                 
           
            
            

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
