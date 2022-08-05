clear
clc
cd /homes/eerfani/riming/SCPP

plot_control
  grey = [0.4,0.4,0.4] ;
  pink = [1.0,0.4,0.6] ;
  purple = [0.5,0,0.5] ;
  phos = [0,0.7,0.7] ;
  violet = [0.2,0.2,0.5] ;
  banafsh = [0.5,0,1] ;
  blu = [0,0.5,1] ;
  grin = [0,1,0.4] ;
      
    fnames = dir('SCPP_all-data_85-87.xlsx') ;

for kk=1:1
  
    [data,txt,raw] = xlsread(fnames(kk).name); 
    LENGTH = data(:,7).*1E+3;  % also converted from mm to um
    mass = data(:,9);  % in mg
    habit = txt(30:end,12) ;

    idx_u = find(ismember(habit,'G1') | ismember(habit,'G3') | ismember(habit,'G3') | ismember(habit,'G5') ...
        | ismember(habit,'G6') ...
        | ismember(habit,'C1C') |ismember(habit,'C1C-A') | ismember(habit,'C1E') ...
        | ismember(habit,'C1E-A') | ismember(habit,'C1F') | ismember(habit,'C1F-A') | ismember(habit,'C2A') ...
        | ismember(habit,'C2B')| ismember(habit,'CP1A') | ismember(habit,'CP1B') | ismember(habit,'CP2A') ...
        | ismember(habit,'CP2B') | ismember(habit,'C2B-A') | ismember(habit,'CP3D') ...
        | ismember(habit,'N1A') | ismember(habit,'N1B') | ismember(habit,'N1C') | ismember(habit,'N1D') ...
        | ismember(habit,'N1E') | ismember(habit,'N1E-A') | ismember(habit,'N2A') | ismember(habit,'N2C') ...
        | ismember(habit,'N2C-A')...
        | ismember(habit,'P1A') | ismember(habit,'P1A-A') | ismember(habit,'P1A-E') | ismember(habit,'P1A-F') ...
        | ismember(habit,'P1B') | ismember(habit,'P1B-F') | ismember(habit,'P1B-FA') | ismember(habit,'P1C') ...
        | ismember(habit,'P1C-A') | ismember(habit,'P1C-F') | ismember(habit,'P1C-FA') | ismember(habit,'P1D') ...
        | ismember(habit,'P1D-F') | ismember(habit,'P1E') | ismember(habit,'P1E-A') | ismember(habit,'P1E-F') ...
        | ismember(habit,'P1E-FA') | ismember(habit,'P1F') | ismember(habit,'P1F-F') | ismember(habit,'P2A') ...
        | ismember(habit,'P2A-A') | ismember(habit,'P2A-F') | ismember(habit,'P2B') | ismember(habit,'P2C') ...
        | ismember(habit,'P2C-F') | ismember(habit,'P2E') | ismember(habit,'P2F') | ismember(habit,'P2F-F') ...
        | ismember(habit,'P2F-A') | ismember(habit,'P2G') | ismember(habit,'P4A') | ismember(habit,'P5') ...
        | ismember(habit,'P7A') | ismember(habit,'P7A-A') | ismember(habit,'P7A-F') | ismember(habit,'P7A-FA') ...
        | ismember(habit,'S1') | ismember(habit,'S1-A') | ismember(habit,'S1-F') | ismember(habit,'S2') ...
        | ismember(habit,'S3') | ismember(habit,'S3-A') ...
        | ismember(habit,'I3A') | ismember(habit,'I3A-F') ) ;
 
      idx_r = find(ismember(habit,'I3B') | ismember(habit,'I3B-F') ...
          | ismember(habit,'R-C1E') | ismember(habit,'R-C1E-A') | ismember(habit,'R-C2B') | ismember(habit,'R-C2B-A') ...
          | ismember(habit,'R-CP1A') ...
          | ismember(habit,'R1D-FA') | ismember(habit,'R2A') | ismember(habit,'R2A-A') | ismember(habit,'R2B') ...
          | ismember(habit,'R2B-F') | ismember(habit,'R2B-FA') | ismember(habit,'R3A') | ismember(habit,'R3B') ...
          | ismember(habit,'R3B-A') | ismember(habit,'R3B-F') | ismember(habit,'R3B-FA') | ismember(habit,'R3C') ...
          | ismember(habit,'R4A') | ismember(habit,'R4B') | ismember(habit,'R4B-A') | ismember(habit,'R-P1A') ...
          | ismember(habit,'R-P1B') | ismember(habit,'R-P1B-F') | ismember(habit,'R-P1B-FA') | ismember(habit,'R-P1C') ...
          | ismember(habit,'R-P1C-F') | ismember(habit,'R-P1C-FA') | ismember(habit,'R-P1D-F') | ismember(habit,'R-P1E') ...
          | ismember(habit,'R-P1E-A') | ismember(habit,'R-P1E-F') | ismember(habit,'R-P1E-FA') | ismember(habit,'R-P2B') ...
          | ismember(habit,'R-P2C') | ismember(habit,'R-P2D') | ismember(habit,'R-P2D-F') | ismember(habit,'R-P2F') ...
          | ismember(habit,'R-P2G') | ismember(habit,'R-P7A') | ismember(habit,'R-P7A-A') | ismember(habit,'R-P7A-FA') ...
          | ismember(habit,'R-N1A') | ismember(habit,'R-N1B') | ismember(habit,'R-N1E') | ismember(habit,'R-N1E-A') ...
          | ismember(habit,'R-N1E-FA') | ismember(habit,'R-N2A') | ismember(habit,'R-N2C') | ismember(habit,'R-N2C-A') ...
          | ismember(habit,'R-S1') | ismember(habit,'R-S1-A') | ismember(habit,'R-S3') ) ;   
    
end


%%%%%%%% plotting    %%%%%%  
 
%%%%% SCPP 
        fig_name = 'bars_dots_m_D_SCPP_r_u';
        fig_dum = figure(1);
      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');


h=plot(LENGTH(idx_u),mass(idx_u),'.');
set(h,'MarkerEdgeColor','b',...
      'MarkerFaceColor','b',...
      'LineWidth',0.6,'MarkerSize',10)
  
hold on
h2=plot(LENGTH(idx_r),mass(idx_r),'.');
set(h2,'MarkerEdgeColor','r',...
      'MarkerFaceColor','r',...
      'LineWidth',0.6,'MarkerSize',10)
  
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Ice Particle Mass (mg)','fontSize',h_axis+6);
box on
ylim([1E-6 1.5E0])
xlim([1E1 1E4])
  

  set(gca,'XMinorTick','on','YMinorTick','on');
  set(gca,'Fontsize',25,'linewidth',1.5)
  set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);

hold on

L_max = 25*1E-4:1E-4:0.6 ; %max(L_tot) ;

% Heymsfield et al.,(2010)
      alpha = 0.00700 ;   % for cgs units (g/cm3)
      beta = 2.20 ; 
      m_h10 = alpha .* L_max .^ beta ;

h3 = loglog(L_max(L_max >= 0.0052 & L_max <=0.3).*1E4,m_h10(L_max >= 0.0052 & L_max <=0.3).*1E3,'-g','LineWidth',3) ;     

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
  

 crvft2 = -0.15980.*log(L_max).^2 + 1.17421.*log(L_max) - 6.72924 ;
h5 = loglog(L_max.*1E4,2.71828.^(crvft2).*1E3,'-k','LineWidth',3) ;
    
hold on
m_sphere = 0.917 .* pi .* L_max .^ 3 ./ 6 ;
m_sphere(m_sphere>0.0005) = NaN ;
h10 = loglog(L_max.*1E4,m_sphere.*1E3,'color',grey,'LineWidth',3) ;

   hleg1 = legend([h h2 h3 h4 h5 h10],'unrimed SCPP data','rimed SCPP data','Heymsfield et al. (2010)',...
       'Cotton et al. (2012)','CPI & cold-habit SCPP curve fit','ice spheres');
    set(hleg1,'Location','SouthEast','Fontsize',19)%h_legend-4)
    set(hleg1,'Interpreter','none')

%equation=sprintf('ice clouds');
cor = '-40\circC < T \leq 0\circC' ;
text = ' ' ;
textt = 'N_r_i_m_e_d_ _&_ _u_n_r_i_m_e_d = 3781' ;
rm = 'N_r_i_m_e_d = 1440' ;
urm = 'N_u_n_r_i_m_e_d = 2341' ;
str = {cor, text, urm, text, rm};
rrr=annotation('textbox', [.15 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',19)
      set(gca,'fontsize',h_axis+6,'LineWidth',2);
      eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
      eval(['print -r600 -dpdf ', fig_name,'.pdf']);       
