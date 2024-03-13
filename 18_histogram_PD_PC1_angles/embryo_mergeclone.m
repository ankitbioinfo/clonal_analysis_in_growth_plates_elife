

clear all 

%path1='2colorCoded3d_for_PC1_and_PD_embryo_uniclone/';
path1='15_micron_2colorCoded3d_for_PC1_and_PD_merged_mask/'

df_file={'S153_m7_distalfemur.dat', 'S154_m3_distalfemur.dat',  'S154_m4_distalfemur.dat'};
pt_file={'S153_m7_proximaltibia.dat', 'S154_m3_proximaltibia.dat',  'S154_m4_proximaltibia.dat'};




df=readfiles(df_file,path1);
pt=readfiles(pt_file,path1);


h=figure;
set(gcf, 'PaperSize', [5 3]); %7
set(gcf, 'PaperPosition', [0 0 5 3]);
bins=30;
histogram(df,bins,'Normalization','pdf')
hold on 
histogram(pt,bins,'Normalization','pdf')
legend('E18.5 DF mergeclone','E18.5 PT mergeclone')
xlabel('[Angle between PD and PC1]')
ylabel('P(PD-PC1)')

saveas(h,['embryo_mergeclone_fill'])
saveas(h,['embryo_mergeclone_fill','.png'])


h=figure;
set(gcf, 'PaperSize', [5 3]); %7
set(gcf, 'PaperPosition', [0 0 5 3]);
histogram(df,bins,'Normalization','pdf','DisplayStyle','stairs' )
hold on 
histogram(pt,bins,'Normalization','pdf','DisplayStyle','stairs')
legend('E18.5 DF mergeclone','E18.5 PT mergeclone')
xlabel('[Angle between PD and PC1]')
ylabel('P(PD-PC1)')

saveas(h,['embryo_mergeclone_staircase'])
saveas(h,['embryo_mergeclone_staircase','.png'])



function d=readfiles(filename,path)
d=[];
for i=1:length(filename)
    %strcat(path,filename{i})
    data=load([path,filename{i}]);
    d=[d;data(:,1)];
end

end