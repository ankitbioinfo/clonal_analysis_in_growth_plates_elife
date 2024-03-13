

clear all 

path1='random_rep2_2colorCoded3d_for_PC1_and_PD/';
path2='convexhull_Embryo_real_2colorCoded3d_for_PC1_and_PD/';

df_file={'S153_m7_distalfemur.dat', 'S154_m3_distalfemur.dat',  'S154_m4_distalfemur.dat'};

%pt_file={'S153_m7_proximaltibia.dat', 'S154_m3_proximaltibia.dat',  'S154_m4_proximaltibia.dat'};




df=readfiles(df_file,path1);
pt=readfiles(df_file,path2);


h=figure;
set(gcf, 'PaperSize', [5 3]); %7
set(gcf, 'PaperPosition', [0 0 5 3]);
bins=30;
histogram(df,bins,'Normalization','pdf')
hold on 
histogram(pt,bins,'Normalization','pdf')
legend('E18.5 DF random','E18.5 DF real')
xlabel('[Angle between PD and PC1]')
ylabel('P(PD-PC1)')

saveas(h,['embryo_rep2_DF_fill'])
saveas(h,['embryo_rep2_DF_fill','.png'])


h=figure;
set(gcf, 'PaperSize', [5 3]); %7
set(gcf, 'PaperPosition', [0 0 5 3]);
histogram(df,bins,'Normalization','pdf','DisplayStyle','stairs' )
hold on 
histogram(pt,bins,'Normalization','pdf','DisplayStyle','stairs')
legend('E18.5 DF random','E18.5 DF real')
xlabel('[Angle between PD and PC1]')
ylabel('P(PD-PC1)')

saveas(h,['embryo_rep2_DF_staircase'])
saveas(h,['embryo_rep2_DF_staircase','.png'])



function d=readfiles(filename,path)
d=[];
for i=1:length(filename)
    %strcat(path,filename{i})
    data=load([path,filename{i}]);
    d=[d;data(:,1)];
end

end