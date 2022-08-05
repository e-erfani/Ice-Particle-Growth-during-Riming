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
  
    
    Min=100;
    Max=4000;
    cint=100;
    bin=[Min:cint:1000 1200 1400 1800  2400 3000 Max];
    bin2=[Min+cint/2:cint:950 1100 1300 1600  2100 2700 3500];

%%%%% SCPP 
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
      
      L_u = LENGTH(idx_u) ; 
      m_u = mass(idx_u) ; 
      L_r = LENGTH(idx_r) ; 
      m_r = mass(idx_r) ; 
    
    for i=1:length(bin)-1
      index=find(L_u<bin(i+1) & L_u>=bin(i));
      if(length(index) < 3)
        mass_new{i} = NaN;
        length_new{i} = NaN;
      else
        mass_new{i}=m_u(index);
        length_new{i}=L_u(index);        
      end
        mass_mean(i)=mean(mass_new{1,i});  
        mass_std(i)=std(mass_new{1,i});
        length_mean(i)=mean(length_new{1,i});
        length_std(i)=std(length_new{1,i});
    end    
    
    for i=1:length(bin)-1
      index=find(L_r<bin(i+1) & L_r>=bin(i));
      if(length(index) < 3)
        mass_new_r{i} = NaN;
        length_new_r{i} = NaN;
      else
        mass_new_r{i}=m_r(index);
        length_new_r{i}=L_r(index);        
      end
        mass_mean_r(i)=mean(mass_new_r{1,i});  
        mass_std_r(i)=std(mass_new_r{1,i});
        length_mean_r(i)=mean(length_new_r{1,i});
        length_std_r(i)=std(length_new_r{1,i});
    end      
end

%%% BL2006 m-A
    fnames2 = dir('BL2006.xlsx') ;

for kk=1:1
    [data2,txt2,raw2] = xlsread(fnames2(kk).name); 
    L_BL = data2(:,4) .* 1E-4;   % converted to cm
    A_BL = data2(:,6) .* 1E-6;  % in mm^2
    m_BL = 0.115 .* A_BL .^ 1.218 ;  % BL2006 m-A: A in mm^2 and m in mg
end

%%%%%%%% plotting    %%%%%%  
        fig_name = 'bars_dots_m_D_SCPP_all';
        fig_dum = figure(1);
      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');

h2=plot(LENGTH(idx_r),mass(idx_r),'.');
set(h2,'MarkerEdgeColor','m',...
      'MarkerFaceColor','m',...
      'LineWidth',0.6,'MarkerSize',9)

hold on  
h=plot(LENGTH(idx_u),mass(idx_u),'.');
set(h,'MarkerEdgeColor','b',...
      'MarkerFaceColor','b',...
      'LineWidth',0.6,'MarkerSize',9)
  
hold on
h3=errorbar(length_mean,mass_mean,mass_std,'ko','LineWidth',4);
set(h3,'MarkerEdgeColor','k',...
      'MarkerFaceColor','k',...
      'LineWidth',1.5)
  
hold on
h4=errorbar_x(length_mean,mass_mean,length_std,'ko');%,'black','LineWidth',3);
set(h4,'MarkerEdgeColor','k',...
      'MarkerFaceColor','k',...
      'LineWidth',1.5)  

hold on
h5=errorbar(length_mean_r,mass_mean_r,mass_std_r,'ro','LineWidth',4);
set(h5,'MarkerEdgeColor','r',...
      'MarkerFaceColor','r',...
      'LineWidth',1.5)
  
hold on
h6=errorbar_x(length_mean_r,mass_mean_r,length_std_r,'ro');%,'black','LineWidth',3);
set(h6,'MarkerEdgeColor','r',...
      'MarkerFaceColor','r',...
      'LineWidth',1.5)  
  
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Ice Particle Mass (mg)','fontSize',h_axis+6);
box on
ylim([1E-4 1.5E0])
xlim([1E2 1E4])
  
hold on
L_max = 150*1E-4:1E-4:0.6 ; %max(L_tot) ;
   
hold on
m_sphere = 0.917 .* pi .* L_max .^ 3 ./ 6 ;
m_sphere(m_sphere>0.001) = NaN ;
h10 = loglog(L_max.*1E4,m_sphere.*1E3,'color',grey,'LineWidth',3) ;
      
hold on
h7=plot(L_BL.*1E4,m_BL,'g*');%,'black','LineWidth',3);
set(h7,'MarkerEdgeColor','g',...
      'MarkerFaceColor','g',...
      'LineWidth',4)

  set(gca,'XMinorTick','on','YMinorTick','on');
  set(gca,'Fontsize',25,'linewidth',1.5)
  set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);

   hleg1 = legend([h h2 h3 h5 h7 h10],'unrimed dendrites','rimed dendrites',...
       'mean unrimed dendrites','mean rimed dendrites','Baker & Lawson 2006','ice spheres');
    set(hleg1,'Location','SouthEast','Fontsize',19)%h_legend-4)
    set(hleg1,'Interpreter','none')

%equation=sprintf('ice clouds');
cor = '-20\circC < T \leq -10\circC' ;
text = ' ' ;
textt = 'N_r_i_m_e_d_ _&_ _u_n_r_i_m_e_d = 970' ;
rm = 'N_r_i_m_e_d = 852' ;
urm = 'N_u_n_r_i_m_e_d = 118' ;
str = {cor, text, urm, text, rm};
rrr=annotation('textbox', [.145 .815, .1, .1], 'String', str);
set(rrr,'Fontsize',18)
      set(gca,'fontsize',h_axis+6,'LineWidth',2);      
errorbar_tick(h3, 0);
herrorbar_tick(h4, 0);
errorbar_tick(h5, 0);
herrorbar_tick(h6, 0);
      
      eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
      eval(['print -r600 -dpdf ' , fig_name,'.pdf']);       
