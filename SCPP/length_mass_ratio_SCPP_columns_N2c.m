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

fnames = dir('SCPP_all-data_85-87.xlsx') ;

for kk=1:1
  
    [data,txt,raw] = xlsread(fnames(kk).name); 
    LENGTH = data(:,7).*1E+3;  % also converted from mm to um
    mass = data(:,9);  % in mg
    habit = txt(30:end,12) ;

    idx_u = find(ismember(habit,'N2C') | ismember(habit,'N2C-A') | ismember(habit,'N2C-F') | ismember(habit,'N2C-FA') ) ;
    
    idx_u_N1e = find(ismember(habit,'N1E') | ismember(habit,'N1E-A') | ismember(habit,'N1E-F') | ismember(habit,'N1E-FA') ) ;
   
      idx_r = find(ismember(habit,'R4C') | ismember(habit,'R4C-A') | ismember(habit,'R4C-F') | ismember(habit,'R4C-FA') ...
          | ismember(habit,'R4B') | ismember(habit,'R4B-A') | ismember(habit,'R4B-F') | ismember(habit,'R4B-FA') ) ;   

      L_u = LENGTH(idx_u) ; 
      m_u = mass(idx_u) ; 
      L_u_N1e = LENGTH(idx_u_N1e) ; 
      m_u_N1e = mass(idx_u_N1e) ; 
      L_r = LENGTH(idx_r) ; 
      m_r = mass(idx_r) ; 
    
    for i=1:length(bin)-1
      index=find(L_u<bin(i+1) & L_u>=bin(i));
      inx_u_1{i} = index ;
        mass_new{i}=m_u(index);
        length_new{i}=L_u(index);        

        mass_mean(i)=mean(mass_new{1,i});  
        mass_std(i)=std(mass_new{1,i});
        length_mean(i)=mean(length_new{1,i});
        length_std(i)=std(length_new{1,i});
    end    

    for i=1:length(bin)-1
      index=find(L_u_N1e<bin(i+1) & L_u_N1e>=bin(i));
      inx_u_1_N1e{i} = index ;
        mass_new_N1e{i}=m_u_N1e(index);
        length_new_N1e{i}=L_u_N1e(index);        

        mass_mean_N1e(i)=mean(mass_new_N1e{1,i});  
        mass_std_N1e(i)=std(mass_new_N1e{1,i});
        length_mean_N1e(i)=mean(length_new_N1e{1,i});
        length_std_N1e(i)=std(length_new_N1e{1,i});
    end    
    
    for i=1:length(bin)-1
      index=find(L_r<bin(i+1) & L_r>=bin(i));
      inx_1{i} = index ;
         mass_new_r{i}=m_r(index);
        length_new_r{i}=L_r(index);        

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

%%%%%% mass ratio plot %%%%

%%%%%%%%%%  R4b & R4c vs. N1e
%%% first, calculate the number of rimed and unrimed particles
for i = 1:length(inx_1)
    n_r(i) = length(inx_1{i}) ;  
end

for i = 1:length(inx_u_1)
    n_u(i) = length(inx_u_1{i}) ;
end
for i = 1:length(inx_u_1_N1e)
    n_u_N1e(i) = length(inx_u_1_N1e{i}) ;
end

%%% Dr. Mitchell's method: weighted average of mass ratio
ratio = mass_mean_r ./ mass_mean ;
ratio_unrimed = mass_mean ./ mass_mean_N1e ;
%mass_mean(mass_mean<0.001) = NaN ;
n_u(ratio>11) = NaN ;
n_u_N1e(ratio>11) = NaN ;
n_r(ratio>11) = NaN ;
ratio(ratio>11) = NaN ;
n_u(n_u==0) = NaN ;
n_u_N1e(n_u_N1e==0) = NaN ;
n_r(n_r==0) = NaN ;

n = 0 ;
for i = 1:length(inx_1)
   if (n_r(i) >= 1)
       n = n + 1 ;
       index_1(n) = i ;
   end   
end

n_u(isnan(n_r)==1) = NaN ;
n_r(isnan(n_u)==1) = NaN ;

ratio_ave = nansum(n_r .* ratio) ./ nansum(n_r) ;

%%% My method: weighted average each rimed and unrimed mass, and then
%%% calculate mass ratio
mass_r_ave = nansum(n_r .* mass_mean_r) ./ nansum(n_r) ;
mass_u_ave = nansum(n_u .* mass_mean) ./ nansum(n_u) ;

ratio_ave_new = mass_r_ave ./ mass_u_ave ;

%%% My method: same as above, but for the shared bins b/w rimed and unrimed
% mass_mean_1_u_r = mass_mean(index_1) ;
mass_r_ave_2 = nansum(n_r .* mass_mean_r) ./ nansum(n_r) ;
mass_u_ave_2 = nansum(n_u(index_1) .* mass_mean(index_1)) ./ nansum(n_u(index_1)) ;

ratio_ave_new_2 = mass_r_ave_2 ./ mass_u_ave_2 ;

%%% My method: {sigma [(nr*mr) / (nu*mu)]} / (Nr / Nu) for the shared bins b/w rimed and unrimed
ratio_ave_new_3 = nansum(n_r .* mass_mean_r ./ (n_u .* mass_mean ) ) ./ (nansum(n_r ./ n_u) ) ; % 
%ratio_ave_new_3 = nansum((n_r .* mass_mean_r ./ (n_u .* mass_mean ) ) ./ (nansum(n_r ./ n_u) ))  ; % 
%ratio_ave_new_3 = nansum(n_r .* mass_mean_r ./ (n_u .* mass_mean ) ) ./ (nansum(n_r) ./ nansum(n_u) ) ; % ./ nansum(n_u_2(index_2)) ) ;


%%%%
  for i = 1:length(bin)*2-2
       bin_PSD(i) = bin(floor(i/2)+1) ;
       ratio_PSD(i) = ratio(round(i/2)) ;       
       ratio_unrimed_PSD(i) = ratio_unrimed(round(i/2)) ;       
  end
    
        fig_name = 'mass_ratio_SCPP_column_N2c';
        fig_dum = figure(2);
      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');
      
h20 = plot(bin_PSD,ratio_PSD,'color',banafsh,'LineWidth',3) ;
hold on
%h30 = plot(bin_PSD,ratio_PSD_2,'--','color',pink,'LineWidth',3) ;

%  set(gca, 'XScale', 'log')
xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('m_r/m_u','fontSize',h_axis+6);
box on
ylim([1 9.5])
%xlim([3E2 1E4])
      set(gca,'XMinorTick','on','YMinorTick','on');
  set(gca,'Fontsize',25,'linewidth',1.5)
  set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);

h100 = scatter(bin2,ratio,'.b','LineWidth',1) ;
for i=1:length(bin2)
    text(bin2(i),ratio(i),num2str(n_r(i)),'horiz','center','vert','bottom','Fontsize',10)
    text(bin2(i),ratio(i),num2str(n_u(i)),'horiz','center','vert','top','Fontsize',10)
end

equation='     -40\circC \leq T \leq -0\circC';
equation2= '' ;
cor=  'Unrimed: N1e' ;
cor2= 'Rimed:    R4b, R4c' ;
str = {equation, equation2, cor, cor2};
rrr=annotation('textbox', [.165 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',18)

      eval(['print -r600 -dpdf ', fig_name,'.pdf']);       
      eval(['print -r600 -djpeg ', fig_name,'.jpg']);   
      
for j = 1:length(bin_PSD)    
    mean_ratio(j) = ratio_ave ;
end
        plot(bin_PSD,mean_ratio,'--','color',pink,'LineWidth',2)

      eval(['print -r600 -dpdf ', fig_name,'ave.pdf']);       
      eval(['print -r600 -djpeg ', fig_name,'ave.jpg']);   

%%%%%%%%
        fig_name = 'u_mass_ratio_SCPP_N2c_N1e';
        fig_dum = figure(10);
      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');
      
h20 = plot(bin_PSD,ratio_unrimed_PSD,'color',banafsh,'LineWidth',3) ;
hold on

xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('m_r/m_u','fontSize',h_axis+6);
box on
ylim([1 2.3])
%xlim([3E2 1E4])
      set(gca,'XMinorTick','on','YMinorTick','on');
  set(gca,'Fontsize',25,'linewidth',1.5)
  set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);

      eval(['print -r600 -dpdf ', fig_name,'.pdf']);       
      eval(['print -r600 -djpeg ', fig_name,'.jpg']);   
