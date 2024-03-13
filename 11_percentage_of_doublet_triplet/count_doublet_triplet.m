
clear all 

%datafile='DF/';
datafile='PT/';
%datafile='P40/';

data=load([datafile,'spherical_coordinate.txt']);

factor=180/pi;
EL=data(:,3)*factor;

count=0;
column_like=[];
for i=1:length(EL)
    if abs(EL(i))>60
        count=count+1;
        column_like(count,:)=data(i,[5,6]);
    end
end

doublet_col=[100*count/length(EL)]

col_unique=length(unique(column_like));
all_unique=length(unique(data(:,[5,6])));

cells_in_doublet=[100*col_unique/all_unique]

LCC=LargestConnectedComponents(column_like);

col2=0;
col3=0;
col4=0;
col5=0;
col6=0;
col7=0;
col8=0;
col9=0;
col10=0;
col11=0;
col12=0;
for i=1:length(LCC)
    if length(LCC{i})==2
        col2=col2+1;
    end
    
    if length(LCC{i})==3
        col3=col3+1;
    end
    
    if length(LCC{i})==4
        col4=col4+1;
    end
    
    if length(LCC{i})>10
        col11=col11+1;
    end
    
     if length(LCC{i})==5
        col5=col5+1;
     end
    
    if length(LCC{i})==6
        col6=col6+1;
    end
    
    if length(LCC{i})==7
        col7=col7+1;
    end
    
    if length(LCC{i})==8
        col8=col8+1;
    end
    
    if length(LCC{i})==9
        col9=col9+1;
    end
    
    if length(LCC{i})==10
        col10=col10+1;
    end
    
end
        
collen=100*[col2,col3,col4,col5,col6,col7,col8,col9,col10,col11]/length(LCC)



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
