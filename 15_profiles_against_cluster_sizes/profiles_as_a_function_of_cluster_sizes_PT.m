clear all 

legname={'E18.5 PT','P40 PT'};


savename={'Volume_normalized','Volume','SurfaceArea_normalized','SurfaceArea','clusterSize','clusterSphericity',  ...
            'clusterPC1','clusterPC2','clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2','radiusOfGyration_normalized','radiusOfGyration',...
            'volumefraction','LorderPar1','LorderPar2','LorderPar3','GorderPar1','GorderPar2','GorderPar3','plane_of_division_growthAxis',...
            'plane_of_division_to_PC1', 'plane_of_division_to_PC2', 'plane_of_division_to_PC3','highest_mode', 'smallest_mode','cluster_mask_fitted_ellipsoid_Rg',...
          'biaxial_S','biaxial_P','biaxial_D','biaxial_C',  'AngleBetweenClusterPC1AndBone_PD' ,'coordination_number','diameter','allometric', 'randomWalktheory','no_ofloops','columnHeight'};

ylabelname={'$\langle\overline{cluster\; volume}\rangle$', '<cluster volume>',  '$\langle\overline{cluster\; surface\; area}\rangle$',  '<cluster surface area>',  '<cluster size>', '<cluster sphericity>',...
  '<cluster PC1>','<cluster PC2>','<cluster PC3>','<cluster PC2/PC1>','<cluster PC3/PC1>','<cluster PC3/PC2>', '$\langle\overline{R_g}\rangle$','<R_g>',...
   '<volume fraction>', '<Local OOP 1>','<Local OOP 2>','Local <OOP 3>',...
   '<Global OOP 1>','<Global OOP 2>','<Global OOP 3>',  '<angle(plane, P-D axis)>','<angle(plane,PC1)>',...
            '<angle(plane, PC2)>', '<angle(plane,PC3)>', '<largest eigenvalue of Hessian>', '<smallest eigenvalue of Hessian>',...
             '<R_g>', '< OOP S>', '< OOP P>', '< OOP D>', '< OOP C>', '<\alpha>','$\langle{avg\;deg}\rangle$','<diameter>','Vol^{2/3}/SA','<R_g>','<# of cycles>','<column height>'};
        
         
NAME={'','just_volume','','just_SA','','sphericity','clusterPC1','clusterPC2','clusterPC3', 'clusterPC2_by_PC1'  ,  'clusterPC3_by_PC1' , 'clusterPC3_by_PC2',...
   '', 'just_rg','VolumeFraction','LOP','LOP','LOP','GOP','GOP','GOP','','','','','','','','','','','',...
    'AngleBetweenClusterPC1AndBone_PD','CluAvgDeg','','allometric','','','colheight'};          
Prop=[2,4,6:12,14:21,33,34,36,39];



plot_clusterProperties=1;

allpathEmbryo={
'data/Nuclei_and_CellsE185_S153_m7_distalfemur/',
'data/Nuclei_and_CellsE185_S153_m7_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m3_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m3_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m4_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m4_proximaltibia/',
}; 

allpathPostnatal={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 
         

dir2=strcat('profiles_against_sizes_PT','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end


embryo_dir='../../EmbryoData/Femur_Tibia/';

path=allpathEmbryo{2};s=strsplit(path,'Nuclei_and_CellsE185_'); p11=strcat('MakeListNucleiLabelled/',s{2});r11=strcat('RandomSampling/',s{2});
path=allpathEmbryo{4};s=strsplit(path,'Nuclei_and_CellsE185_'); p12=strcat('MakeListNucleiLabelled/',s{2});r12=strcat('RandomSampling/',s{2});
path=allpathEmbryo{6};s=strsplit(path,'Nuclei_and_CellsE185_'); p13=strcat('MakeListNucleiLabelled/',s{2});r13=strcat('RandomSampling/',s{2});
rad1=load([embryo_dir,p11,'ParameterForRandomModel_basedon_nuc_.dat']);
rad2=load([embryo_dir,p12,'ParameterForRandomModel_basedon_nuc_.dat']);
rad3=load([embryo_dir,p13,'ParameterForRandomModel_basedon_nuc_.dat']);
embryoradius={rad1,rad2,rad3};


p40_dir='../proximalTibia_DistalFemur/';
path=allpathPostnatal{2};s=strsplit(path,'Nuclei_and_CellsP40_'); pn11=strcat('MakeListNucleiLabelled/',s{2});r11=strcat('RandomSampling/',s{2});
path=allpathPostnatal{4};s=strsplit(path,'Nuclei_and_CellsP40_'); pn12=strcat('MakeListNucleiLabelled/',s{2});r12=strcat('RandomSampling/',s{2});
path=allpathPostnatal{6};s=strsplit(path,'Nuclei_and_CellsP40_'); pn13=strcat('MakeListNucleiLabelled/',s{2});r13=strcat('RandomSampling/',s{2});
rad1=load([p40_dir,pn11,'ParameterForRandomModel_basedon_nuc_.dat']);
rad2=load([p40_dir,pn12,'ParameterForRandomModel_basedon_nuc_.dat']);
rad3=load([p40_dir,pn13,'ParameterForRandomModel_basedon_nuc_.dat']);
p40radius={rad1,rad2,rad3};


if plot_clusterProperties==1
        

% Property is function of cluster size 
% h=figure;
% XL=0.06;XR=0.02;XGap=0.06;Row=4;
% YT=0.06;YB=0.1;YGap=0.04;Col=5;
% Width=(1-XL-XR-((Col-1)*XGap))/Col;
% Height=(1-YT-YB-((Row-1)*YGap))/Row;
% YPos=1-YT-Height; 
% set(gcf, 'PaperSize', [15 7]); %7
% set(gcf, 'PaperPosition', [0 0 15 7]);

    



for chro=1:length(NAME)
   
        h=figure;
        set(gcf, 'PaperSize', [5 3]); %7
        set(gcf, 'PaperPosition', [0 0 5 3]);
            
               PZzone=2;  %it is between 0.1 to 0.8
               embryo_data_path={p11,p12,p13};
               %rand_data_path={r11,r12,r13};
               data=[]; Xbin=[];     
               for real_samp=1:length(embryo_data_path)
                    d3=load([embryo_dir, embryo_data_path{real_samp},'AllFeaturesSave_.mat']);
                       if Prop(chro)==16
                           value=d3.(NAME{Prop(chro)})(:,1);
                       elseif Prop(chro)==17
                           value=d3.(NAME{Prop(chro)})(:,2);
                       elseif Prop(chro)==18
                           value=d3.(NAME{Prop(chro)})(:,3);
                       elseif Prop(chro)==19
                           value=d3.(NAME{Prop(chro)})(:,1);
                       elseif Prop(chro)==20
                           value=d3.(NAME{Prop(chro)})(:,2);
                       elseif Prop(chro)==21
                           value=d3.(NAME{Prop(chro)})(:,3);
                       else
                           value=d3.(NAME{Prop(chro)});
                       end    
                        data=[data; d3.clusterSize real(value)]; 
                        Xbin=[Xbin; embryoradius{real_samp}(1:length(d3.clusterSize),1)];
               end
               [p(1),ln,ylimit(1,:)]=FunctionOfCellSize(data, Xbin,2,3,11,0,PZzone,'r.-'); 
               % title('cre clusters','fontweight','normal') 
              
              
                p40_data_path={pn11,pn12,pn13};
                 data=[]; Xbin=[];     
                for real_samp=1:length(p40_data_path)
                    d3=load([p40_dir, p40_data_path{real_samp},'AllFeaturesSave_.mat']);
                       if Prop(chro)==16
                           value=d3.(NAME{Prop(chro)})(:,1);
                       elseif Prop(chro)==17
                           value=d3.(NAME{Prop(chro)})(:,2);
                       elseif Prop(chro)==18
                           value=d3.(NAME{Prop(chro)})(:,3);
                       elseif Prop(chro)==19
                           value=d3.(NAME{Prop(chro)})(:,1);
                       elseif Prop(chro)==20
                           value=d3.(NAME{Prop(chro)})(:,2);
                       elseif Prop(chro)==21
                           value=d3.(NAME{Prop(chro)})(:,3);
                       else
                           value=d3.(NAME{Prop(chro)});
                       end    
                        data=[data; d3.clusterSize real(value)]; 
                        Xbin=[Xbin; p40radius{real_samp}(1:length(d3.clusterSize),1)];
               end
               [p(2),ln,ylimit(2,:)]=FunctionOfCellSize(data, Xbin,2,3,11,0,PZzone,'b.-'); 
             
              
      
%                data=[]; Xbin=[];
%                for rand_samp=1:length(rand_data_path)
%                for realization=1:no_of_realization
%                     d3=load([rand_data_path{rand_samp},'AllFeaturesSave',num2str(realization),'.mat']);
%                     
%                       if Prop(chro)==16
%                            value=d3.(NAME{Prop(chro)})(:,1);
%                        elseif Prop(chro)==17
%                            value=d3.(NAME{Prop(chro)})(:,2);
%                        elseif Prop(chro)==18
%                            value=d3.(NAME{Prop(chro)})(:,3);
%                        elseif Prop(chro)==19
%                            value=d3.(NAME{Prop(chro)})(:,1);
%                        elseif Prop(chro)==20
%                            value=d3.(NAME{Prop(chro)})(:,2);
%                        elseif Prop(chro)==21
%                            value=d3.(NAME{Prop(chro)})(:,3);
%                        else
%                            value=d3.(NAME{Prop(chro)});
%                       end    
%                     
%                     
%                     data=[data; [d3.clusterSize real(value)]];  
%                     Xbin=[Xbin; radius{rand_samp}(1:length(d3.clusterSize),1)];
%                end
%                end
%                
%             
%                [p(2),ln,ylimit(2,:)]=FunctionOfCellSize(data, Xbin,2,3,11,0,PZzone,'b.-'); 

            %FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,Prop(chro),1,2);
            
            ylim([min(ylimit(:,1)), max(ylimit(:,2))])
            if chro==19
                   ylabel(ylabelname{Prop(chro)},'fontsize',8,'Interpreter','Latex');
            else
                   ylabel(ylabelname{Prop(chro)},'fontsize',8);
            end
            
                
                legend(p,legname,'location','best','fontsize',5);
                
            

            
                xlabel('cluster size (# of cells)')
                xlim([1,50]);

    %chro            
    %NAME{Prop(chro)}
    %savename{chro}
    saveas(h,[dir2,savename{Prop(chro)}])
    saveas(h,[dir2,savename{Prop(chro)},'.png'])

    close all 

end

end