clear all 

allpath={
'data/Nuclei_and_CellsE185_S153_m7_distalfemur/',
'data/Nuclei_and_CellsE185_S153_m7_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m3_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m3_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m4_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m4_proximaltibia/',
};     


Area=[];
for gi=[2,4,6]%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        inputpath=strcat('MakeListNucleiLabelled/',s{2});
        data=load([inputpath,'AllFeaturesSave_.mat']);
        %outputpath=strcat('RealSampling/',s{2});
        Area=[Area;data.just_volume];
end     






g = fittype( @(k, x) (k^k)*( x.^(k-1)).* exp(-k.*x) / gamma(k));

scale1 = Area - min(Area);
scale2 = mean(Area) -min(Area);
ratio=scale1/scale2;
r=ratio(find(ratio<3));
%r=ratio;

h=figure;
set(gcf, 'PaperSize', [4 3]); %7
set(gcf, 'PaperPosition', [0 0 4 3]);
[yr,xr]=histnorm(r);
bar(xr,yr);
hold on 

[fE,GE] = fit(xr',yr',g,'StartPoint',[1]);
%p8=plot(xr,fE(xr),'k--','Linewidth',2);

coef=fE.k;
%coef=0.1;
plot(xr,g(coef,xr),'k--','Linewidth',2);


xlabel('$\frac{Vol - min(Vol)}  {\overline{Vol} - min(Vol)}$' , 'Interpreter','Latex');
ylabel('P()')
		
title(['distal femur: fit k = ',num2str(coef)])
saveas(h,'gamma_.png')
%print(gcf,'gamma_.png')
