clear
clc
plot_control
  grey = [0.4,0.4,0.4] ;
  pink = [1.0,0.4,0.6] ;
  purple = [0.5,0,0.5] ;
  phos = [0,0.7,0.7] ;
  violet = [0.2,0.2,0.5] ;
  banafsh = [0.5,0,1] ;
  blu = [0,0.5,1] ;
  grin = [0,1,0.4] ;
  

for kk=1:1
  
    [data,txt,raw] = xlsread(fnames(kk).name); 
    LENGTH = data(:,7).*1E+3;  % also converted from mm to um
    mass = data(:,9);  % in mg

    Min=100;
    Max=4000;
    cint=100;
    bin=[Min:cint:1000 1200 1400 1800  2400 3000 Max];
    bin2=[Min+cint/2:cint:950 1100 1300 1600  2100 2700 3500]; 
    
    for i=1:length(bin)-1
      index=find(LENGTH<bin(i+1) & LENGTH>=bin(i));        
        mass_new{i}=mass(index);
        mass_mean(i)=mean(mass_new{1,i});  
        mass_std(i)=std(mass_new{1,i});
        length_new{i}=LENGTH(index);
        length_mean(i)=mean(length_new{1,i});
        length_std(i)=std(length_new{1,i});
    end    

end


%%%%%%%% plotting    %%%%%%  
    fnames = dir('SCPP_all-data_85-87.xlsx') ;

        fig_name = 'bars_fit_m_D_SCPP_all';
        fig_dum = figure(1);
      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');
    

L_max = 25*1E-4:1E-4:0.4 ; %max(L_tot) ;

% Heymsfield et al.,(2010)
      alpha = 0.00700 ;   % for cgs units (g/cm3)
      beta = 2.20 ; 
      m_h10 = alpha .* L_max .^ beta ;

h3 = loglog(L_max(L_max >= 0.0052 & L_max <=0.3).*1E4,m_h10(L_max >= 0.0052 & L_max <=0.3).*1E3,'g','LineWidth',3) ;     

  hold on

% Cotton et al.,(2012)
      alph1 = (pi/6.)*0.700 ;   % for cgs units (g/cm3)
      beta1 = 3.00 ;
      alph2 = 0.00257 ; % for cgs units
      beta2 = 2.00 ;
      
      m_c12(L_max <= 0.007) = alph1 .* L_max(L_max <= 0.007) .^ beta1 ;
      m_c12(L_max > 0.007) = alph2 .* L_max(L_max > 0.007) .^ beta2 ;

h4 = loglog(L_max(L_max >= 0.003 & L_max <=0.05).*1E4,m_c12(L_max >= 0.003 & L_max <=0.05).*1E3,'--k','LineWidth',3) ;    
  hold on
  
  
%%%%% SCPP 
h=errorbar(length_mean,mass_mean,mass_std,'r.');%,'black','LineWidth',3);
set(h,'MarkerEdgeColor',banafsh,...
      'MarkerFaceColor',banafsh,...
      'LineWidth',1.5)
  
hold on
h2=errorbar_x(length_mean,mass_mean,length_std,'r.');%,'black','LineWidth',3);
set(h2,'MarkerEdgeColor',banafsh,...
      'MarkerFaceColor',banafsh,...
      'LineWidth',1.5)
  
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Ice Particle Mass (mg)','fontSize',h_axis+6);
box on
ylim([1E-6 1E0])
xlim([1E1 1E4])
  
errorbar_tick(h, 0);
herrorbar_tick(h2, 0);


  set(gca,'XMinorTick','on','YMinorTick','on');
  set(gca,'Fontsize',25,'linewidth',1.5)
  set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);

hold on
 crvft2 = -0.15980.*log(L_max).^2 + 1.17421.*log(L_max) - 6.72924 ;
h5 = loglog(L_max.*1E4,2.71828.^(crvft2).*1E3,'-k','LineWidth',3) ;
    
hold on
m_sphere = 0.917 .* pi .* L_max .^ 3 ./ 6 ;
m_sphere(m_sphere>0.0005) = NaN ;
h10 = loglog(L_max.*1E4,m_sphere.*1E3,'color',grey,'LineWidth',3) ;

   hleg1 = legend([h h3 h4 h5 h10],'mean SCPP data','Heymsfield et al. (2010)',...
       'Cotton et al. (2012)','CPI & cold-habit SCPP curve fit','ice spheres');
    set(hleg1,'Location','SouthEast','Fontsize',19)%h_legend-4)
    set(hleg1,'Interpreter','none')

%equation=sprintf('ice clouds');
cor = '-40\circC < T \leq 0\circC' ;
text = ' ' ;
textt = 'N_S_C_P_P = 4869' ;
str = {cor , text, textt};
rrr=annotation('textbox', [.15 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',19)
      set(gca,'fontsize',h_axis+6,'LineWidth',2);
      eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
      eval(['print -r600 -dpdf ', fig_name,'.pdf']);       

      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Percent Difference
 
curve = exp(-0.15980.*log(length_mean.*1E-4).^2 + 1.17421.*log(length_mean.*1E-4) - 6.72924) .* 1E3 ;
 diff = 100 .* (curve - mass_mean) ./ ((curve + mass_mean) ./ 2) ;
 
         fig_name = 'percent_diff';
        fig_dum = figure(2);
      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');
      
      h100 = scatter(length_mean,diff)   ;
   set(h100,'MarkerEdgeColor',purple,...
      'MarkerFaceColor',purple,...
      'LineWidth',4)
 
  hold on
  L_max = 25*1E-4:1E-4:1 ;
  zero = 0 .* L_max ;
  plot(L_max.*1E4,zero,'k--')
  
set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
xlabel('Ice Particle Mean Size (\mum)','fontSize',h_axis+6);
ylabel('Percent Difference','fontSize',h_axis+6);
box on
ylim([-40 40])
xlim([1E2 1E4])
  set(gca,'XMinorTick','on','YMinorTick','on');
  set(gca,'Fontsize',25,'linewidth',1.5)
  set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);

  cor = '-40\circC < T \leq 0\circC' ;
text = ' ' ;
textt = 'Comparing mean masses using all SCPP particles' ;
textt2 = 'with corresponding masses from curve fit' ;
str = {textt, textt2};
%rrr=annotation('textbox', [.15 .2, .1, .1], 'String', str);
%set(rrr,'Fontsize',19)
  
  set(gca,'fontsize',h_axis+6,'LineWidth',2);
      eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
