clear all 
 
 
no_of_realization=20;

plot_rg=0;
plot_cono=0;
plot_fractal=1;
plot_clusterProperties=0;



savename={'Volume_normalized','Volume','SurfaceArea_normalized','SurfaceArea','clusterSize','clusterSphericity',  ...
            'clusterPC1','clusterPC2','clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2','radiusOfGyration_normalized','radiusOfGyration',...
            'volumefraction','LorderPar1','LorderPar2','LorderPar3','GorderPar1','GorderPar2','GorderPar3','plane_of_division_growthAxis',...
            'plane_of_division_to_PC1', 'plane_of_division_to_PC2', 'plane_of_division_to_PC3','highest_mode', 'smallest_mode','cluster_mask_fitted_ellipsoid_Rg',...
          'biaxial_S','biaxial_P','biaxial_D','biaxial_C',  'AngleBetweenClusterPC1AndBone_PD' ,'coordination_number','diameter','allometric', 'randomWalktheory','no_ofloops'};

ylabelname={'$\langle\overline{cluster\; volume}\rangle$', '<cluster volume>',  '$\langle\overline{cluster\; surface\; area}\rangle$',  '<cluster surface area>',  '<cluster size>', '<cluster sphericity>',...
  '<cluster PC1>','<cluster PC2>','<cluster PC3>','<cluster PC2/PC1>','<cluster PC3/PC1>','<cluster PC3/PC2>', '$\langle\overline{R_g}\rangle$','<R_g>',...
   '<volume fraction>', '<Local OOP 1>','<Local OOP 2>','Local <OOP 3>',...
   '<Global OOP 1>','<Global OOP 2>','<Global OOP 3>',  '<angle(plane, P-D axis)>','<angle(plane,PC1)>',...
            '<angle(plane, PC2)>', '<angle(plane,PC3)>', '<largest eigenvalue of Hessian>', '<smallest eigenvalue of Hessian>',...
             '<R_g>', '< OOP S>', '< OOP P>', '< OOP D>', '< OOP C>', '<\alpha>','$\langle{avg\;deg}\rangle$','<diameter>','Vol^{2/3}/SA','<R_g>','<# of cycles>'};
        
         
latexSymbol=[1,3,13,34];         
        

dir2=strcat('RealRandomCluster_profiles','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end


allpath={
'data/Nuclei_and_CellsE185_S153_m7_distalfemur/',
'data/Nuclei_and_CellsE185_S153_m7_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m3_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m3_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m4_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m4_proximaltibia/',
};         
        

count=1;
for fi=[1,3,5]
    path=allpath{fi};s=strsplit(path,'Nuclei_and_CellsE185_'); 
    dataFile1=strcat('MakeListNucleiLabelled/',s{2});
    real_data{count}=load([dataFile1,'SaveVariablesForPlot_.mat']);
    
    dataFile1=strcat('RandomSampling/',s{2});
    random_data{count}=load([dataFile1,'SaveVariablesForPlot.mat']);
    count=count+1;
end
         
for k=1:length(savename)
      clear temp 
      for i=1:3
           df=real_data{i};
           temp(:,i)=df.AvgMeanFeatures{k};           
      end
      
      mu=nanmean(temp,2);
      sd=nanstd(temp')';
      realFeat{k}=[mu,sd];
      X=df.Xbin; 
end


%          
%          for i=1:3
%      df=d_femur{i};
%      dt=d_tibia{i};
%      df_data(:,i)=df.Ybin;
%      dt_data(:,i)=dt.Ybin;
% end
% 
%  v1=nanmean(df_data,2);  sd1=nanstd(df_data')';
%  v2=nanmean(dt_data,2);  sd2=nanstd(dt_data')';
% 
% 
% dataFile={'Tibia/','NullModel/'};
% mycolor={'ro-','bs-'};
% 
% d={real_data,random_data};


n=size(real_data{1}.dataAvgMeanFeatures,1);
m=size(real_data{1}.dataAvgMeanFeatures,2);

for i=1:n % 38 properites 
    for j=1:m % xbin data 
            data1=real_data{1}.dataAvgMeanFeatures{i,j};   % individual point 
            data2=(random_data{1}.realization_AvgMeanFeatures{i}(j,:))';
            [size(data1), size(data2)]
            [h,p]=ttest2(data1,data2);
            pvaluetest{i}(j,1)=p;
            
    end
end

radius=load('Tibia/dataSave/ParameterForRandomModel_basedon_nuc_medium.dat');
radius=load([outputpath_labelled,'ParameterForRandomModel_basedon_nuc_.dat']);



if plot_cono==1

% coordination number 
h=figure;
XL=0.09;XR=0.02;XGap=0.06;Row=2;
YT=0.06;YB=0.1;YGap=0.09;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 
set(gcf, 'PaperSize', [6 4]); %7
set(gcf, 'PaperPosition', [0 0 6 4]);
zones={'RZ','PZ','PHZ','HZ'};
for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        subplot('Position',marray);
        
        d3=load([dataFile{1},'dataSave/AllFeaturesSave_medium.mat']);
        data=[d3.clusterSize d3.coordNumber]; Xbin=radius(1:length(d3.clusterSize),1);
        [p(1),ln]=FunctionOfCellSize(data, Xbin,2,3,31,1,chro,'r.'); lname{1}=['cre : ', ln];
        %FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,31,1,chro,'r.');
        hold on 
        
        data=[]; Xbin=[]; 
        for realization=1:no_of_realization
            d3=load([dataFile{2},'Realization_medium/AllFeaturesSave',num2str(realization),'.mat']);
            data=[data; [d3.clusterSize d3.coordNumber]];
            Xbin=[Xbin; radius(1:length(d3.clusterSize),1)];
        end
                
        [p(2),ln]=FunctionOfCellSize(data, Xbin,2,3,31,1,chro,'b.'); lname{2}=['background : ', ln];
        
        legend(p, lname,'location','northwest','fontsize',5)
        
        title(zones{chro},'fontweight','normal')
        if j==1
        ylabel(ylabelname{31},'fontsize',9)
        end
        if i==Row
        xlabel('cluster size (# of cells)')
        else
            set(gca,'xticklabel',[])
        end
        ylim([1.2,7])
        xlim([0,40])
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end
%saveas(h,[dir2,'CoordinationNumber.png'])
print([dir2,'comparison_CoordinationNumber.png'],'-dpng','-r300')
close all 
end




if plot_rg==1
%Radius of gyration from Random walk model 
h=figure;
XL=0.09;XR=0.02;XGap=0.06;Row=2;
YT=0.06;YB=0.1;YGap=0.09;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 
set(gcf, 'PaperSize', [6 4]); %7
set(gcf, 'PaperPosition', [0 0 6 4]);
zones={'RZ','PZ','PHZ','HZ'};
for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        subplot('Position',marray);
        
        
        d3=load([dataFile{1},'dataSave/AllFeaturesSave_medium.mat']);
        data=[d3.clusterSize d3.rg]; Xbin=radius(1:length(d3.clusterSize),1);
        [p(1),ln]=FunctionOfCellSize(data, Xbin,2,3,11,0,chro,'r.-'); 
        hold on 
        data=[d3.clusterSize d3.averageStepSize];
        [p(2),ln]=FunctionOfCellSize(data, Xbin,2,3,33,0,chro,'k-'); 

        
        
        data=[]; Xbin=[]; data2=[];
        for realization=1:no_of_realization
            d3=load([dataFile{2},'Realization_medium/AllFeaturesSave',num2str(realization),'.mat']);
            data=[data; [d3.clusterSize d3.rg]];
            data2=[data2; [d3.clusterSize d3.averageStepSize]];
            Xbin=[Xbin; radius(1:length(d3.clusterSize),1)];
        end
                
        [p(3),ln]=FunctionOfCellSize(data, Xbin,2,3,11,0,chro,'b.-'); 
        [p(4),ln]=FunctionOfCellSize(data2, Xbin,2,3,33,0,chro,'k--');
        
       % legend(p, lname,'location','northwest','fontsize',5)
       lname={'cre cluster', 'RW theory for cre','background cluster','RW theory for background'};
        

%         p(1)=FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,11,0,chro);
%         hold on 
%         p(2)=FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,33,0,chro);
        if chro==1
        legend(p,lname,'location','northwest','fontsize',5)
        end
        title(zones{chro},'fontweight','normal')
        if j==1
        ylabel(ylabelname{11})
        end
        if i==Row
        xlabel('cluster size (# of cells)')
        else
            set(gca,'xticklabel',[])
        end
        ylim([4.5,28])
        xlim([0,40])
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end
print([dir2,'comparison_RadiusOfGyration.png'],'-dpng','-r300')
close all 
end






% fractal like pattern 
if plot_fractal==1
 h=figure;
XL=0.09;XR=0.02;XGap=0.05;Row=1;
YT=0.11;YB=0.15;YGap=0.12;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 
set(gcf, 'PaperSize', [7 3]); %7
set(gcf, 'PaperPosition', [0 0 7 3]);
for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        subplot('Position',marray);
        
        if chro==1
               d3=load([dataFile{1},'dataSave/AllFeaturesSave_medium.mat']);
               
               %average step size is basically theoretical RW so real step
               %size l is defined as   l = sqrt(6/n) Rg 
               l=sqrt(6./d3.clusterSize).*d3.averageStepSize;                
               data=[d3.clusterSize d3.rg]; Xbin=radius(1:length(d3.clusterSize),1);
                fractalCalculation( data, Xbin,l)
                title('cre clusters: R_g = pf*l*N^{\alpha(=1/D_f)}','fontweight','normal') 
                ylabel('log(R_g)'); 
        end
        
        if chro==2
               data=[]; Xbin=[];lstep=[];
               for realization=1:no_of_realization
                   d3=load([dataFile{2},'Realization_medium/AllFeaturesSave',num2str(realization),'.mat']);
                   l=sqrt(6./d3.clusterSize).*d3.averageStepSize; 
                   data=[data; [d3.clusterSize d3.rg]];  
                   Xbin=[Xbin; radius(1:length(d3.clusterSize),1)];
                   lstep=[lstep; l];
               end
                fractalCalculation( data, Xbin,lstep)
                title('background clusters: R_g = pf*l*N^{\alpha(=1/D_f)}','fontweight','normal') 
        end
        
        ylim([1.3,3.3])
        xlim([1,4])
        
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end
%saveas(h,[dir2,'fractalDimension.png'])
print([dir2,'comparison_fractalDimension.png'],'-dpng','-r300')
close all 
end



pause 




 
 

 
 



h=figure;
set(gcf, 'PaperSize', [5 3]); %7
set(gcf, 'PaperPosition', [0 0 5 3]);
for k=1:length(dataFile)
    if k==1
    q(k)=plot(d1.Xbin,d1.Ybin,mycolor{k},'linewidth',1);  
    hold on    
    else
         q(k)=errorbar(d2.Xbin,d2.Ybin,d2.YbinStd,mycolor{k},'linewidth',1); 
    end
    
        
    
end

legname={'cre cluster','background cluster'};

legend(q,legname,'location','northeast')
xlabel('Long axis of bone')
ylabel('cluster density');
saveas(h,[dir2,'clusterDensity.png'])
close all 




for k=1:length(savename)
    h=figure;
    set(gcf, 'PaperSize', [5 3]); %7
    set(gcf, 'PaperPosition', [0 0 5 3]);
    for i=1:length(dataFile)
             data=d{i};
             index=1:length(data.Xbin);                   
             q(i)=errorbar(data.Xbin(index),data.AvgMeanFeatures{k}(index),data.AvgStdFeatures{k}(index),mycolor{i},'linewidth',1);
            
            hold on 
    end
    
   %legend(q,legname,'location','northeast')  %'northeast'
    
    
                   for j=1:length(pvaluetest{k})
                       if pvaluetest{k}(j,1)<0.05
                           plot(d{1}.Xbin(j),d{1}.AvgMeanFeatures{k}(j),'ks','linewidth',2,'markersize',15,'markerfacecolor','none');
                       end
                   end
    
               
    
    if (k>=13)&(k<=15)
        title('local','fontweight','normal')
        legend(q,legname,'location','south')  %'northeast'
        ylim([0,1.01])
    end
    
     if (k>=16)&(k<=18)
        title('global','fontweight','normal')
        legend(q,legname,'location','south')
        ylim([-0.69,1])
     end
     
     if (k>=19)&(k<=22)
	if k==22
        	legend(q,legname,'location','northwest')  %'northeast'
 	else
		legend(q,legname,'location','northeast')  %'northeast'
	end
         ylim([10,80])
     end
    
    
    
    
                  
    
   
    
    xlim([0,1])
    xlabel('Long axis of bone')
    ylabel(ylabelname{k});
    saveas(h,[dir2,savename{k},'.png'])
    close all 
end






if plot_clusterProperties==1

    ylabelname={'<volume>', '<surface area>','<size>', '<sphericity>',...
  '<PC1 coeff>','<PC2 coeff>','<PC3 coeff>','<PC2/PC1 coeff>','<PC3/PC1 coeff>','<PC3/PC2 coeff>', '<R_g>',...
   '<\phi>', '<Local OOP 1>','<Local OOP 2>','Local <OOP 3>',...
   '<Global OOP 1>','<Global OOP 2>','<Global OOP 3>',  '<angle(plane, P-D axis)>','<angle(plane,PC1)>',...
            '<angle(plane, PC2)>', '<angle(plane,PC3)>', '<largest eigenvalue of Hessian>', '<smallest eigenvalue of Hessian>',...
             '<cluster radius>', '< OOP S>', '< OOP P>', '< OOP D>', '< OOP C>', '<\alpha>','<coordination number>','<diameter>'};

% Property is function of cluster size 
h=figure;
XL=0.06;XR=0.02;XGap=0.06;Row=4;
YT=0.06;YB=0.1;YGap=0.04;Col=5;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 
set(gcf, 'PaperSize', [15 7]); %7
set(gcf, 'PaperPosition', [0 0 15 7]);
Prop=[1,2,4,5:18,30,31,32];
NAME={'volume','surface_area','','sphericity','clusterPC1','clusterPC2','clusterPC3', 'clusterPC2_by_PC1'  ,  'clusterPC3_by_PC1' , 'clusterPC3_by_PC2',...
    'rg','VolumeFraction','LOP','LOP','LOP','GOP','GOP','GOP','','','','','','','','','','','',...
    'AngleBetweenClusterPC1AndBone_PD','coordNumber','diameter'}; 
    



for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        
        if chro<=length(Prop)
            subplot('Position',marray);
            
               PZzone=2;
               d3=load([dataFile{1},'dataSave/AllFeaturesSave_medium.mat']);
               
               if Prop(chro)==12
                   value=d3.(NAME{Prop(chro)})(:,1);
               elseif Prop(chro)==13
                   value=d3.(NAME{Prop(chro)})(:,2);
               elseif Prop(chro)==14
                   value=d3.(NAME{Prop(chro)})(:,3);
               elseif Prop(chro)==15
                   value=d3.(NAME{Prop(chro)})(:,1);
               elseif Prop(chro)==16
                   value=d3.(NAME{Prop(chro)})(:,2);
               elseif Prop(chro)==17
                   value=d3.(NAME{Prop(chro)})(:,3);
               else
                   value=d3.(NAME{Prop(chro)});
               end    
                   
      
               data=[d3.clusterSize real(value)]; Xbin=radius(1:length(d3.clusterSize),1);
               [p(1),ln,ylimit(1,:)]=FunctionOfCellSize(data, Xbin,2,3,11,0,PZzone,'r.-'); 
               % title('cre clusters','fontweight','normal') 
              
               
      
               data=[]; Xbin=[];
               for realization=1:no_of_realization
                    d3=load([dataFile{2},'Realization_medium/AllFeaturesSave',num2str(realization),'.mat']);
                    
                      if Prop(chro)==12
                           value=d3.(NAME{Prop(chro)})(:,1);
                       elseif Prop(chro)==13
                           value=d3.(NAME{Prop(chro)})(:,2);
                       elseif Prop(chro)==14
                           value=d3.(NAME{Prop(chro)})(:,3);
                       elseif Prop(chro)==15
                           value=d3.(NAME{Prop(chro)})(:,1);
                       elseif Prop(chro)==16
                           value=d3.(NAME{Prop(chro)})(:,2);
                       elseif Prop(chro)==17
                           value=d3.(NAME{Prop(chro)})(:,3);
                       else
                           value=d3.(NAME{Prop(chro)});
                      end    
                    
                    
                    data=[data; [d3.clusterSize real(value)]];  
                    Xbin=[Xbin; radius(1:length(d3.clusterSize),1)];
               end
               
            
               [p(2),ln,ylimit(2,:)]=FunctionOfCellSize(data, Xbin,2,3,11,0,PZzone,'b.-'); 

            %FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,Prop(chro),1,2);
            
            ylim([min(ylimit(:,1)), max(ylimit(:,2))])
            ylabel(ylabelname{Prop(chro)},'fontsize',8)
            
            lname={'cre cluster','background cluster'};
            if chro==1
                legend(p,lname,'location','northwest','fontsize',5);
            end

            if i==Row 
                xlabel('cluster size (# of cells)')
            else
                set(gca,'xticklabel',[])
            end

        end
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end
saveas(h,[dir2,'VariousMetricFunctionOfClusterSize.png'])
close all 








% coordination number 
% real_size= mydata{3,j};
% real_rg=mydata{11,j};
% avg_deg=mydata{31,j};
% phi=mydata{12,j};
%  h=figure; 
% CoordinatinNumberCalculation(d1.dataAvgMeanFeatures, d1.Xbin,1,12,31)
% xlabel('Volume fraction')
% ylabel('Coordination number');
% saveas(h,[dir2,'phi_vs_avgCoordinatinNumber.png'])
% close all 
% 
% 
%  
% 
% 
% 
% 
%  
%  h=figure;
%  plot(d1.Xbin, d1.AvgMeanFeatures{31}, 'r.-');
%  hold on 
%  errorbar(d1.Xbin,d1.AvgMeanFeatures{31},d1.AvgStdFeatures{31},'r','linewidth',1);
%  ylabel('Coordination number');     
%  xlabel('Long axis of bone')
%  saveas(h,[dir2,'SpatialProfile_AvgCoordinatinNumber.png'])
%  close all 

end






function fractalCalculation( mydata, Xbin,stepsize)  
        
    g=fittype(@(m,c,x) (m*x+c));
    RZ=[];
    PZ=[];
    PHZ=[];
    HZ=[];
    for j=1:length(Xbin)
       
       real_size=mydata(j,1);
       real_rg=mydata(j,2);
       RW_l=stepsize(j,1);
        
            
        if Xbin(j)<=0.2
            RZ=[RZ; [log(real_size),log(real_rg),RW_l]];
        elseif (Xbin(j)>0.2)&(Xbin(j)<=0.5)
            PZ=[PZ; [log(real_size),log(real_rg),RW_l]];
        elseif (Xbin(j)>0.5)&(Xbin(j)<=0.7)
            PHZ=[PHZ; [log(real_size),log(real_rg),RW_l]];
        elseif (Xbin(j)>0.7)&(Xbin(j)<=1)
            HZ=[HZ; [log(real_size),log(real_rg),RW_l]];
        end 

    end
    
    zone={RZ,PZ,PHZ,HZ};
    zonecolor={'r.','b.','g.','k.'};
    legname={'RZ','PZ','PHZ','HZ'};
    
    for i=1:length(zone)
        data=zone{i};
        t=min(data(:,1)): 0.2: max(data(:,1))+0.000001;
        %t=unique(data(:,1));
        
        X=[];
        Y=[];
        Z=[];
        L=[];
        for j=1:length(t)-1
            index=find((data(:,1)>= t(j) ) & (  data(:,1)<t(j+1)));
            %index=find(data(:,1)==t(j));
            X=[X; t(j)];
            Y=[Y;  nanmean(data(index,2))];
            Z=[Z;  nanstd(data(index,2))];
            L=[L;  nanmean(data(index,3))];
        end            
         
        index = find(~isnan(Y));
        x=X(index); y=Y(index); 
        plot(x,y,zonecolor{i});
        hold on 
        
        %[x,y]
        
        min_y1=y-Z(index); max_y1=y+Z(index);  index=1:length(x);
        stacky2=(min_y1);stacky1=(max_y1);
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
        h=fill(fillxx,fillyy,zonecolor{i}(1),'EdgeColor','none','facealpha',0.2);
        
        
       
        [f,~]=fit(x,y,g,'startpoint',[1,0]);
        p(i)=plot(x,f(x),[zonecolor{i}(1),'-'],'linewidth',0.5);

        Df=1/f.m; 
        slope=f.m;
        
        % R_g = pf*n^alpha*l
        % log(R_g) = alpha*log(N) + log(l) + log(pf) 
        
        pf =  exp(f.c)/nanmean(L);
        
        lname{i}=[legname{i},': ', '\alpha =',sprintf('%0.2f',slope), ', D_f =',sprintf('%0.2f',Df), ', pf=', sprintf('%0.2f',pf)];
    
    end

        legend(p, lname,'location','southeast','fontsize',5)
        xlabel('log(N)');
        

   
end



function Area=CalculateSurfaceArea(K,v)
          Area=0;
          for surfind=1:size(K,1)
              pointA=v(K(surfind,1),:);
              pointB=v(K(surfind,2),:);
              pointC=v(K(surfind,3),:);
              TriangleVertex=[pointA; pointB; pointC];
              Ts=pdist(TriangleVertex,'euclidean');
              p=sum(Ts)/2;
              Area=Area+sqrt(p*(p-Ts(1))*(p-Ts(2))*(p-Ts(3)));
          end
          
end


function graphEnergy=calculateEnergy(node,cent)
    
    A=zeros(length(node));
    D=zeros(length(node));
    centroid=cent(node,:);
    for i=1:size(centroid,1)
        for j=i+1:size(centroid,1)
            dist=pdist(centroid([i,j],:));
            A(i,j)= 1/dist;
            A(j,i)= 1/dist;
        end
    end


    for i=1:size(A,1)
        D(i,i)=sum(A(i,:));
    end

    normalizedAdjacency=D^(-1/2)*A*(D^(-1/2));
    L=D-A;
    normalizedLaplacian=D^(-1/2)*L*(D^(-1/2));


    E1=eig(normalizedAdjacency);
    E2=eig(normalizedLaplacian);

    %E2([1:3,end])

    energyAdjacency=sum(abs(E1));

    %totalWeight=sum(sum(A));
    totalWeight=sum(E2);

    energyLaplacian=sum( abs( E2 - (totalWeight/size(A,1))));

    secondSmallestLaplacianEigenvalue=E2(2);

    graphEnergy=[ secondSmallestLaplacianEigenvalue,energyAdjacency,energyLaplacian,min(E1),max(E1)];


end



