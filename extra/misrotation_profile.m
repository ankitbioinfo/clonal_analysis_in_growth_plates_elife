clear all 


warning('off','all')


allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 


RZ_HZ_height={ [-353,511], [-218,635], [-384,1154],[-232,587],[-445,1004],[-227,473]};


for gi=6%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        outputpath=strcat('MakeListNucleiLabelled/',s{2});
        outputpath_mask=strcat('MakeListClustersMask/',s{2});
        %savedata=strcat('Misrotation/',s{2});
end        

load([outputpath,'centroid_and_surface_nuclei.mat']);
sarah_mask=load([outputpath_mask,'centroid_and_surface_nuclei.mat']);
sarah_mask_id=load([outputpath,'ClusterMask_good.dat']);

[~,LCC]=readClusterFile([outputpath,'Cluster_good.dat']); 


bone_mu=mean(centroid);
c=centroid-bone_mu;

bininterval=[0,45,60,70,92];
mycolorAll={'b','r','g','y'};
mycolorAll={'y','g','r','b'};

misrotationData=[];

h=figure;
fname=s{2}(1:strlength(s{2})-1);
BoneGrowthAxis=[0; 0 ; 1];

fid_col=fopen(['Misrotation/',fname,'_col.txt'],'w');
fid_clu=fopen(['Misrotation/',fname,'_clu.txt'],'w');


%for j=1:5%:length(LCC)
for j=1:length(LCC)
     ClusterShape=[];
     id=LCC{j};
     mergeVol=0;
     for i=1:length(LCC{j})
            %ClusterShape=[  ClusterShape;nuc{LCC{j}(i)}];
            mergeVol=mergeVol+celvolume(LCC{j}(i));
     end
   
     combinedCellVolumeInCluster(j,1)=mergeVol;
     %[K,V]=convhull(ClusterShape); 
     
    mskc=sarah_mask.nuc{sarah_mask_id(j)}; 
    ClusterShape=mskc;
    CC=mean(ClusterShape); 
    [clusterPCvector,~,latent]=pca(mskc);
%     plot3(mskc(:,1),mskc(:,2),mskc(:,3),'k.','markersize',0.1);
%     hold on 

    AngleBetweenClusterPC1AndBone_PD=90-180/pi*oangle(clusterPCvector(:,1),BoneGrowthAxis);
    
    %[min(mskc(:,3)), max(mskc(:,3))]-bone_mu(3)
    
    
    [CluAvgDeg(j,1), ColAvgDeg(j,1),largestDistance(j,1),stepsize,deg,edges,AZ_EL_R,misrot,zaxis]=coordination_number(c,LCC{j}); 
    EL=abs(AZ_EL_R(:,3))*180/pi;
    ELsort=sort(EL,'descend');
    p10=ceil(0.1*length(ELsort));
    p20=ceil(0.2*length(ELsort));
    p50=ceil(0.5*length(ELsort));
    anglechosen=ELsort(p50);
    
    length1= max(mskc(:,3))-min(mskc(:,3));
    length2= max(zaxis(:,3))-min(zaxis(:,3));

    
    for k=1:length(misrot)
        %column id, column length,cellid, AngleBetweenClusterPC1AndBone_PD,  misrot, zaxis 
       
        if AngleBetweenClusterPC1AndBone_PD>60
           fprintf(fid_col,'%d\t%d\t%d\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n'  , j,length(id), k, AngleBetweenClusterPC1AndBone_PD, 180/pi*misrot(k),zaxis(k,3),length1,length2);
        else
           fprintf(fid_clu,'%d\t%d\t%d\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n'  , j,length(id), k, AngleBetweenClusterPC1AndBone_PD, 180/pi*misrot(k),zaxis(k,3),length1,length2);
        end

    end
   
    
    for k=1:length(bininterval)-1
         if (anglechosen>=bininterval(k))&(anglechosen<bininterval(k+1))
                mycolor=mycolorAll{k};
         end
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
% 
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

fclose(fid_clu);
fclose(fid_col);



% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% 
% savepath='./';

 % saveas(h,[savepath,'colorcoded',fname,'.png']);
%  saveas(h,[savepath,'PC1_colormap_50p_',fname,'.fig']);


%dlmwrite('angledist2.dat',[ankit50, ankit60, ankit70])



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


function [avgdeg, avgdegExpectedColumn,largestDistance,averageStepSize,deg,E,sph,misrotation,zaxis]=coordination_number(centroid,nodeid)
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
             
            
            misrot=zeros(length(unique(E(:))),1);
            normalize=zeros(length(unique(E(:))),1);
            zaxis=zeros(length(unique(E(:))),3);    
            normzaxis=zeros(length(unique(E(:))),3); 
            
            count=1;
            for i=1:size(E,1)
                c=cent(E(i,:),:);
                c=c-mean(c);
                [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
                flag=[1;1];
                sph(count:count+1,:)=[flag,az,el,r,nodeid(E(i,[1,2]))',nodeid(E(i,[2,1]))'];
                count=count+2;
               
                misrot(  E(i,1))=misrot(  E(i,1))+abs(el(1));
                misrot(  E(i,2))=misrot(  E(i,2))+abs(el(1));
                normalize(E(i,1))=normalize(  E(i,1))+1;
                normalize(E(i,2))=normalize(  E(i,2))+1;
                zaxis(E(i,2),:)=cent(E(i,2),:);
                zaxis(E(i,1),:)=cent(E(i,1),:);
          
            end    
            
            misrotation=misrot./normalize;
            %[misrot, normalize,misrotation]
           
            
            
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
