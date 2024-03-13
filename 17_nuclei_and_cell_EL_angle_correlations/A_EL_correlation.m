clear all 

dir='Nuclei_and_Cells_DT_S17_m2_wt/';

a=load([dir,'all_cells.mat']);
b=load([dir,'all_cells_nuclei.mat']);

cells=a.all_cells(:,5:7);
nuc=b.all_cells_nuclei(:,5:7);

cpc1=a.all_cells(:,8:10);
cpc2=a.all_cells(:,11:13);
cpc3=a.all_cells(:,14:16);

npc1=b.all_cells_nuclei(:,8:10);
npc2=b.all_cells_nuclei(:,11:13);
npc3=b.all_cells_nuclei(:,14:16);


% each row of all_cells_nuclei is nucleus
% the columns contains the individual nuclei features
% 1 - stack id
% 2 - volume
% 3 - surface area
% 4 - sphericity
% 5-7 - centroid x,y,z coordinates
% 8-10 - PC1 x,y,z orientation
% 11-13 - PC2 x,y,z orientation
% 14-16 - PC3 x,y,z orientation
% 17-19 - PC1,PC2,PC3 latent coefficient
% 20 - Delaunay density

% 
% figure 
% plot3(nuc(:,1),nuc(:,2),nuc(:,3),'b.')
fid=fopen('matching_pairs.dat','w');

for i=1:length(cells)
    for j=1:length(nuc)
        dist=pdist([cells(i,:);nuc(j,:)]);
        if dist<20
            ang1= oangle( cpc1(i,:)', npc1(j,:)');
            ang2= oangle( cpc2(i,:)', npc2(j,:)');
            ang3= oangle( cpc3(i,:)', npc3(j,:)');
            fprintf(fid,'%d\t%d\t%0.3f\t%0.2f\t%0.2f\t%0.2f\n',i,j,dist,ang1,ang2,ang3);
        end
    end
end
fclose(fid);

% 
% badcell=[];
%  for j=1:size(edges,1)
%      
%      dist=pdist(centroid(edges(j,:),:));
%      if dist<8
%              cell1=surfaces{edges(j,1)}; actualVol(1)=volume( edges(j,1)   );      
%              cell2=surfaces{edges(j,2)}; actualVol(2)=volume( edges(j,2)   );  
%              combined=[cell1;cell2];
%              [~,combV]=convhull(combined);
% 
%              volumeFraction= combV/sum(actualVol);
%              if volumeFraction<1
%                         [sa,sb]=min(actualVol);
%                         if sb==1
%                             badid=edges(j,1);
%                         end
%                         if sb==2
%                             badid=edges(j,2);
%                         end
%                         badcell=[badcell;badid];
%              end
%      end
%  end    

%  goodcellindex=setdiff(firstgood,badcell);
%  toc
  %length(goodcellindex)  
  
 function  angle= oangle(u,v)
      angle=atan2(norm(cross(u,v)),dot(u,v));
      if angle>(pi/2)
          angle=pi-angle;
      end
      angle=angle*180/pi;
 end
