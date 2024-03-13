function [p,lname,ylimit]=FunctionOfCellSize( mydata, Xbin,flag,P1,P2,fitflag,myzone,zonecolor)


%    flag=1 for real and flag=2 for random 
%    P1 is x axis and P2 is yaxis features 
%    fitflag whether you want fit or not 
%    myzone which bone zone want to plot 

%    data1=d1.dataAvgMeanFeatures{i,j};   % individual point 
%    data2=(d2.realization_AvgMeanFeatures{i}(j,:))';
%    
%     
        
    
    RZ=[];
    PZ=[];
    PHZ=[];
    HZ=[];
    for j=1:length(Xbin)
        
%         if flag==1
%             real_size= mydata{P1,j};
%             real_rg=mydata{P2,j};
%           
%         end    
        
        if flag==2
            real_size=mydata(j,1);
            real_rg=mydata(j,2);
        end
            
        if Xbin(j)<=0.1
            RZ=[RZ; [real_size,real_rg]];
        elseif (Xbin(j)>0.1)&(Xbin(j)<=0.8)
            PZ=[PZ; [real_size,real_rg]];
        elseif (Xbin(j)>0.8)&(Xbin(j)<=0.85)
            PHZ=[PHZ; [real_size,real_rg]];
        elseif (Xbin(j)>0.85)&(Xbin(j)<=1)
            HZ=[HZ; [real_size,real_rg]];
        end 

    end
    
    
    zone={RZ,PZ,PHZ,HZ};
    %[length(RZ), length(PZ), length(PHZ), length(HZ)]
    %zonecolor={'r.','b.','g.','k.'};
    legname={'RZ','PZ','PHZ','HZ'};
    
    D1=[];
    D2=[];
    
    g=fittype(@(m,c,x) (m*x+c));
    g1=fittype(@(alpha,beta,x) alpha*(1-exp(beta*x)));
    g2=fittype(@(alpha,beta,x) alpha*(exp(beta*x)));
    
    for i=myzone%1:length(zone)
        data=zone{i};
        X=unique(data(:,1));
        clear Y 
        for j=1:length(X)
            yindex=find(data(:,1)==X(j));
            Y(j,1)= mean(data(yindex,2));
            Y(j,2)= std(data(yindex,2));
        end
        
        
        
        
        if P2==33
            
            p=plot(X,Y(:,1),zonecolor,'markersize',5,'linewidth',1,'markerfacecolor',zonecolor(1));
        else
            p=plot(X,Y(:,1),zonecolor,'markersize',5,'linewidth',1,'markerfacecolor',zonecolor(1));
        end
            
        
        hold on 
        %errorbar(X,Y(:,1),Y(:,2),zonecolor{i},'linewidth',0.5);

        d1= Y(:,1)-Y(:,2);  D1=[D1;d1];
        d2= Y(:,1)+Y(:,2);  D2=[D2;d2];
        index=1:(length(X)-1);
        min_y1=d1; max_y1=d2;  index=1:length(X);
        stacky2=(min_y1);stacky1=(max_y1);
        fillxx=X(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
        
        if fitflag==0
            h=fill(fillxx,fillyy,zonecolor(1),'EdgeColor','none','facealpha',0.2);
        else
            h=fill(fillxx,fillyy,zonecolor(1),'EdgeColor','none','facealpha',0.2);
        end
%         ax = gca;
%         ax.YAxis.Exponent = 2;
        lname=0;
        if fitflag==1
                if ((P2==31))
                    [f,~]=fit(X(index),Y(index,1),g1,'startpoint',[1,-1]);
                    Df1=f.alpha; Df2=f.beta;
                    lname=strcat('Y = ' ,sprintf('%0.1f',Df1),'*(1-exp(',sprintf('%0.2f',Df2),'x))');   
                elseif (P2==12)    
                    [f,~]=fit(X(index),Y(index,1),g2,'startpoint',[1,-1]);
                    Df1=f.alpha; Df2=f.beta;
                    lname=strcat('Y = ' ,sprintf('%0.2f',Df1),'*exp(',sprintf('%0.2f',Df2),')');  
                else
                    [f,~]=fit(X(index),Y(index,1),g,'startpoint',[1,0]);
                    Df1=f.m; Df2=f.c;
                    lname=strcat('Y = ',sprintf('%0.2f',Df1), '*x + ',sprintf('%0.2f',Df2)    );
                end
                 p=plot(X,f(X),strcat(zonecolor(1),'-'),'linewidth',1);
                 %title(lname{count},'fontweight','normal','fontsize',8)                     
         end

         
        
    end

   
    %ylim([min(D1), max(D2)]);
    ylimit=[min(D1),max(D2)];
    
end
