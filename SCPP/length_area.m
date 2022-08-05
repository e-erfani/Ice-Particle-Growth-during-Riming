clear
clc

fname='ice_crystals.xls';
c=xlsread(fname);

LENGTH=c(:,5);
area=c(:,7);


figure(1); 
%plot(ave_variable,ave_height,'magenta','LineWidth',3);
h=scatter(LENGTH,area);%,'black','LineWidth',3);
set(h,'MarkerEdgeColor','b',...
      'MarkerFaceColor',[0 .5 .5],...
      'LineWidth',0.6)
%axis([0 360 0 5000])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('ice particle size (microns)','fontSize',26);
ylabel('ice particle area (square microns)','fontSize',26);
%title('Vertical Profile of Wind Direction','fontsize',30,'fontweight','bold');
box on

hold on
p = polyfit(log10(LENGTH),log10(area),1);
f = polyval(p,log10(LENGTH));
loglog(LENGTH,10.^f,'-k')
[correlation,P_value]=corrcoef(log10(LENGTH),log10(area));

cor=sprintf('R=%g',correlation(2,1));
equation=sprintf('\lambda y=%gx^{%g}',10.^p(2),p(1));
str = {equation, cor};
rrr=annotation('textbox', [.2 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',21)
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)


%%

% Min=100;
% Max=4700;
% cint=100;
% bin=Min:cint:Max;
% bin2=Min+cint/2:cint:Max-cint/2;
% 
% for i=1:length(bin)-1
%     index=find(LENGTH<=bin(i+1) & LENGTH>bin(i));        
%     
%     area_new{i}=area(index);
%     area_mean(i)=nanmean(area_new{1,i});
%     area_std(i)=nanstd(area_new{1,i});
%     
% end

a0=-2.20919228;
a1=1.39612329;
a2=-0.0523313433;

L=150:3500;
l_cm=L.*10^(-4);
log_l=log(l_cm);
log_A=a0+a1.*log_l+a2.*log_l.^2;
A_cm=exp(log_A);
A=A_cm.*10.^8;

%ezplot('10^8.*exp(-2.2092+1.3961.*log(x.*10^(-4))-0.0523.*((log(x.*10^(-4))).^2))',[100 10000])

A_sphere=0.25.*pi.*L.^2;

Min=100;
Max=4000;
cint=100;
bin=[Min:cint:1000 1200 1400 1800  2400 3000 Max];
bin2=[Min+cint/2:cint:950 1100 1300 1600  2100 2700 3500]; 

for i=1:length(bin)-1
    index=find(LENGTH<=bin(i+1) & LENGTH>bin(i));        
    
    area_new{i}=area(index);
    area_mean(i)=mean(area_new{1,i});  
%    area_mean(i)=nanmean(area_new{1,i});  
    area_std(i)=std(area_new{1,i});
%    area_std(i)=nanstd(area_new{1,i});
    length_new{i}=LENGTH(index);
    length_mean(i)=mean(length_new{1,i});
%    length_mean(i)=nanmean(length_new{1,i});
    length_std(i)=std(length_new{1,i});
%    length_std(i)=nanstd(length_new{1,i});    
end



figure(2); 
%plot(ave_variable,ave_height,'magenta','LineWidth',3);
h=errorbar(bin2,area_mean,area_std,'r.');%,'black','LineWidth',3);
set(h,'MarkerEdgeColor','b',...
      'MarkerFaceColor',[0 .5 .5],...
      'LineWidth',0.6)
  
hold on
%h2=herrorbar(bin2,area_mean,length_std,'o');%,'black','LineWidth',3);
h2=errorbar_x(bin2,area_mean,length_std,'r.');%,'black','LineWidth',3);
set(h2,'MarkerEdgeColor','b',...
      'MarkerFaceColor',[0 .5 .5],...
      'LineWidth',0.6)
  
%axis([0 360 0 5000])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Ice Particle Size (µm)','fontSize',26);
ylabel('Ice Particle Projected Area (µm^2)','fontSize',26);
%title('Vertical Profile of Wind Direction','fontsize',30,'fontweight','bold');
box on        
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)        
errorbar_tick(h, 0);
herrorbar_tick(h2, 0);

 hold on
 h3=plot(L,A,'-k','linewidth',1.5);
 
 hold on
 h4=plot(L,A_sphere,'--k','linewidth',1.5);

 hleg1 = legend([h h3 h4],'SCPP Data','SPARTICUS Curve fit','Ice Spheres');
    set(hleg1,'Location','NorthWest','Fontsize',20)
    set(hleg1,'Interpreter','none')
    
    
percent_std=100.*area_std./area_mean;
final_perc_std=nanmean(percent_std);    
