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


for gi=2:5%length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        inputpath=strcat('MakeListNuclei_All/',s{2});
        outputpath=inputpath;
       

load([inputpath,'centroid_and_surface_nuclei.mat']);


%[~,LCC]=readClusterFile([inputpath,'Cluster_NoDuplication.dat']); 


[All_edges,min_neigh_dist,neighlist,edges,remaining_edges,cel_normalizationFactor] = calculate_nuclei_density(centroid, [1, 1, 1], 2);
 
dlmwrite([outputpath,'neighbors_edges.dat'],All_edges,'\t');

end



function [All_edges,neighbor,neighborList,edges1,edges2,convexVolume] = calculate_nuclei_density(N, spacing, delta_spacing, exclude_boundary_points)
%{
Input:
    N - an n*3 matrix with the [x,y,z] coordinates of each of the n points.
    spacing -  1*3 array with the [x,y,z] physical size of the voxels (if N is already in physical units, just use [1,1,1]).
    delta_spacing - a positive real value that specifies the frequency of the sampling of the 3D grid in the output image (the density space / the output argument V).
    exclude_boundary_points - a boolean value that allows the user to choose whether to include or exclude boundary vertices to avoid potential noise in the
    physical boundaries of the sample, default is 'false' (i.e. include boundaries)
Output:
    V - the 3D matrix of estimated densities.
    X, Y, Z - the coordinates of the interpolated space V.
Example:
    N = rand(100, 3);
    spacing = [1 1 1];
    delta_spacing = 0.01;
    exclude_boundary_points = false;
    V = calculate_nuclei_density(N, spacing, delta_spacing, exclude_boundary_points);
    figure;
    imagesc(V(:,:,50));
    axis image;
    colormap jet;
%}
% if the user didn't define 'exclude_boundary_points' we set it to false:
if ~exist('exclude_boundary_points', 'var')
    exclude_boundary_points = false;
end
% calculating the triangulation and the volume of each triangle:
TRI = delaunay(N(:,1), N(:,2), N(:,3));
[~,convexVolume]=convexHull(delaunayTriangulation(N));

clear neighbor
All_edges=[];
for i = 1 : size(N,1)
    temp=[];
    for j=1:size(TRI,1)
        for k=1:size(TRI,2)
            if TRI(j,k)==i
                temp=[temp,TRI(j,:)];
            end
        end
    end

   
    neighborList{i}=setdiff(unique(temp),i);
    %neighbor(i,1)=length(neighborList{i});
    ids= neighborList{i};
    
    clear dist 
    for k=1:length(ids)
        dist(k)=pdist(N([i,ids(k)],:)); 
        temp=[i, ids(k), dist(k)];
        All_edges=[All_edges; temp];
    end
   
    [sa,sb]=sort(dist);
    neighbor(i,1)=min(dist);
    
    edges1(i,:) = [sort([i,ids(sb(1))]) sa(1)  ];  
    edges2{i}=[ids(sb(2:end))',sa(2:end)'];
end




%triplot(tri,x,y); for 2d 
%tetramesh(tri,X); 3d 

[size(N,1), length(unique(TRI(:)))];

end
         