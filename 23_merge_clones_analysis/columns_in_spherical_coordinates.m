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

mkdir Fig 

for gi=3%length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        inputpath=strcat('MakeListNuclei_All/',s{2});
        
        data=load([inputpath,'centroid_and_surface_nuclei.mat']);
        
        neighbors=load([inputpath,'neighbors_edges.dat']);
        E=neighbors(:,[1,2]);
        E=sort(E')';
        [~,firstgood]=unique(E,'rows');
        edges=E(firstgood,:);
        dist=neighbors(firstgood,3);
        
        length(edges)
        index=find(dist<17);
        ed=edges(index,:);
        
        length(ed)
        
        sph=coordination_number(data.centroid,ed);
        EL=abs(sph(:,2)*180/pi); 
        index=find(EL>75); 
        sp=sph(index,:);
        %dlmwrite([inputpath,'ElevationAngleSphericalCoordinates.dat'],sp,'\t');
        
        column_like_ed=sp(:,[4,5]); %E=sort(coled')';
        LCC=LargestConnectedComponents(column_like_ed);
        fid=fopen([inputpath,'ExpectColumnLikeCluster.dat'],'w');
        
        for i=1:length(LCC)
            for j=1:length(LCC{i})
                fprintf(fid,'%d ',LCC{i}(j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        
        
end

name=s{2}(1:strlength(s{2})-1);
for i=1:20
    h=figure;
myid=LCC{i};
for fi=1:length(myid)
    c=data.nuc{myid(fi)};
    plot3(c(:,1),c(:,2),c(:,3),'r.','markersize',1)
    hold on 
    cent=data.centroid(myid(fi),:);
    text(cent(:,1),cent(:,2),cent(:,3),num2str(myid(fi)))
end

saveas(h,['Fig/Largest_expected_extracted_columns_',name,num2str(i),'.png']);
saveas(h,['Fig/Largest_expected_extracted_columns_',name,num2str(i),'.fig']);
close all 
end





function [sph]=coordination_number(cent,neighbors)
E=neighbors(:,[1,2]);
 count=1;
 for i=1:size(E,1)
     c=cent(E(i,:),:);
     c=c-mean(c);
     [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
     sph(count:count+1,:)=[az,el,r,(E(i,[1,2]))',(E(i,[2,1]))'];
     count=count+2;
 end
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