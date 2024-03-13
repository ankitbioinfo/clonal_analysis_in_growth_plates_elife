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


RZ_HZ_height={ [-353,511], [-218,635], [-384,1154],[-232,587],[-445,1004],[-227,473]};



for gi=6%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        %outputpath=strcat(randpath,s{2});
        %outputpath_mask=strcat('MakeListClustersMask/',s{2});
        outputpath_unlabelled=strcat('../Femur_Tibia/MakeListNuclei_unlabelled/',s{2});
        rand_dir_output=strcat('../Femur_Tibia/RandomSampling/',s{2});
end        

load([outputpath_unlabelled,'centroid_and_surface_nuclei.mat']);
%sarah_mask=load([outputpath_mask,'centroid_and_surface_nuclei.mat']);
%sarah_mask_id=load([outputpath,'ClusterMask_NoDuplication.dat']);


realization=1;
[LCC,allid]=readClusterFile(strcat(rand_dir_output,'Realization/Random_cluster_',num2str(realization),'.dat'));


bone_mu=mean(centroid);
c=centroid-bone_mu;

bininterval=[0,40,50,60,70,80,92];
%bininterval=[0,10,20,30,40,50,92];
%mycolorAll={'b','r','m','g','c','y'};
mycolorAll={'y','c','g','m','r','b'};
BoneGrowthAxis=[0; 0 ; 1];



section=unique(unique_tileid(:,3));

for sec=1:length(section) 
    %for j=1:5%:length(LCC)
    h=figure;

for j=1:length(LCC)
     ClusterShape=[];
     id=LCC{j};
     mergeVol=0;
     for i=1:length(LCC{j})
            ClusterShape=[  ClusterShape;nuc{LCC{j}(i)}];
            mergeVol=mergeVol+celvolume(LCC{j}(i));
     end
   
     combinedCellVolumeInCluster(j,1)=mergeVol;
     %[K,V]=convhull(ClusterShape); 
     
    %mskc=sarah_mask.nuc{sarah_mask_id(j)}; 
    %ClusterShape=mskc;
    CC=mean(ClusterShape); 
    
    [clusterPCvector,~,latent]=pca(ClusterShape);
    
%     plot3(mskc(:,1),mskc(:,2),mskc(:,3),'k.','markersize',0.1);
%     hold on 
    
    
   % [CluAvgDeg(j,1), ColAvgDeg(j,1),largestDistance(j,1),stepsize,deg,edges,AZ_EL_R]=coordination_number(c,LCC{j}); 
         
    AngleBetweenClusterPC1AndBone_PD=90-180/pi*oangle(clusterPCvector(:,1),BoneGrowthAxis);
   
    data(j,1)=AngleBetweenClusterPC1AndBone_PD;
    data(j,2)=section(sec);
    data(j,3)=unique_tileid(id(1),1);
    
    
    for k=1:length(bininterval)-1
         if (AngleBetweenClusterPC1AndBone_PD>=bininterval(k))&(AngleBetweenClusterPC1AndBone_PD<bininterval(k+1))
                mycolor=mycolorAll{k};
         end
    end
    
      for k=1:length(id)
         ct=unique_tileid(id(k),3);
    end
    
    
    
%     A70=sum(EL>70)/length(EL); 
%     A60=sum(EL>60)/length(EL); 
%     A45=sum(EL>45)/length(EL);
%     if A70 > 0.15
%         mycolor=mycolorAll{1};
%     elseif A60 > 0.3 
%         mycolor= mycolorAll{2};
%     elseif A45 > 0.5
%         mycolor=mycolorAll{3};
%     else
%         mycolor=mycolorAll{4};
%     end

    
        if ct==section(sec)


    plot3(CC(:,1),CC(:,2),CC(:,3),'k.','markersize',5);
    hold on 
    [vec,val]=eig(cov(ClusterShape));
    ovec=vec;
    oval=val;
    mu=mean(ClusterShape);
    d = sqrt(diag(val));
    hold on;
    %factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10); %PC3
    %factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'g','LineWidth',5); %PC2
    factor=1; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),mycolor,'LineWidth',2); %PC1

    %factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10); %PC3 
    %factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'g','LineWidth',5); %PC2 
   factor=-1; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),mycolor,'LineWidth',2); %PC1 
        end
% 

%    veccolor=0;
%     countfreq=zeros(1,length(bininterval)-1);
%     for i=1:length(EL)
%         for k=1:length(bininterval)-1
%             if (EL(i)>=bininterval(k))&(EL(i)<bininterval(k+1))
%                 %veccolor=mycolor{k};
%                 countfreq(k)=countfreq(k)+1;
%             end
%         end
%     end
%     
%     
%     ankit(j,:)=countfreq/length(EL);
%     
%     ankit70(j,:)=sum(EL>70)/length(EL); 
%     ankit60(j,:)=sum(EL>60)/length(EL); 
%     ankit50(j,:)=sum(EL>45)/length(EL); 
    %[j, countfreq, length(EL)]
    

end


xlabel('X')
ylabel('Y')
zlabel('Bone long axis (RZ-HZ)')
title('E18.5')

savepath='./';
fname=s{2}(1:strlength(s{2})-1);

%set(gca,'Color','k');
 % saveas(h,[savepath,'colorcoded',fname,'.png']);
  saveas(h,[savepath,'5colorCoded3d_for_PC1_and_PD/sec_',num2str(section(sec)),'_colormap_between_PD_and_PC1_',fname,'.fig']);
  dlmwrite([savepath,'5colorCoded3d_for_PC1_and_PD/',fname,'.dat'],data);
  close all 
end

colperc=100*sum(data(:,1)>=60)/length(data)
cluperc=100*sum(data(:,1)<60)/length(data)


% scheme 1 
%    if the cluster has > 15% doublets angle cutoff > 70 then it is blue 
%elseif the cluster has > 30% doublets angle cutoff > 60 then it is red 
%elseif the cluster has > 50% doublets angle cutoff > 45 then it is green 
%or 
%elseif the cluster has > 45% doublets angle cutoff > 50 then it is green 


% scheme 2 
%    if the cluster has > 10% doublets angle cutoff > 75 then it is blue 
%elseif the cluster has > 30% doublets angle cutoff > 60 then it is red
%elseif the cluster has > 50% doublets angle cutoff > 45 then it is green 
%othwise yellow 



 function  angle= oangle(u,v)
      angle=atan2(norm(cross(u,v)),dot(u,v));
      if angle>(pi/2)
          angle=pi-angle;
      end
 end


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
            
            %[length(E),length(myedge)]

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