% Author: Ehsan Erfani

clc; clear; close all;
plot_control
cd /homes/eerfani/riming

  grey = [0.4,0.4,0.4] ;
  pink = [1.0,0.4,0.6] ;
  purple = [0.5,0,0.5] ;
  phos = [0,0.7,0.7] ;
  violet = [0.2,0.2,0.5] ;
  banafsh = [0.5,0,1] ;
  blu = [0,0.5,1] ;
  grin = [0,1,0.4] ;

%%%%%%%%%%%%  constants
gg = 981 ; 
rhow = 1.0 ;
RA = 286.7 ; % SI or mks
TK = -8 + 273.15 ; %  temperature in Kelvin units
P = 800 ; % hPa
rhoa=1E-3.*108.7.*P ./ (RA.*TK) ; % g/cm**3
vis= 2.48E-6 .* TK .^ 0.754 ;
visk = vis./rhoa ;
      
%%%%%%%%%   ice particle dimension and m-D & A-D coefficients 
   D_w = 47E-4 : 1E-4 : 292.8E-4 ;
   D_l = 412.81 .* D_w .^2 - 4.5773 .* D_w + 0.0221 ;
   D_r = sqrt(D_l .^ 2 + D_w .^ 2) ;    

     beta(D_r > 0.03) = 1.74 ;  % hex col.
     alpha(D_r > 0.03) = 0.000907 ;
     delta(D_r > 0.03) = 1.414 ;
     gamma(D_r > 0.03) = 0.0512 ;  

     beta(D_r > 0.01 & D_r <= 0.03) = 1.91 ;  % hex col.
     alpha(D_r > 0.01 & D_r <= 0.03) = 0.00166 ;
     delta(D_r > 0.01 & D_r <= 0.03) = 1.50 ;
     gamma(D_r > 0.01 & D_r <= 0.03) = 0.0696 ;  
     
     beta(D_r <= .01) = 2.91 ;  % hex col.
     alpha(D_r <= 0.01) = 0.1677 ;
     delta(D_r <= 0.01) = 2.00 ;
     gamma(D_r <= 0.01) = 0.684 ;       
   
      del=5.83 ; % for ice particles
      c0_0=0.6 ; % for ice particles
      t1i = del.^2 ./ 4 ;
      t2i = (1./t1i) .* (1./c0_0.^0.5) ;
           
%    Define correction terms for turbulence and aggregation:
      a1=0.0017 ;  % applied to aggregates
      b1= 0.8   ;  % applied to aggregates
         
%%%%%%%%%%%%  Best number and Reynols number         
       xx = (2.*alpha.*gg.*rhoa.*D_r.^(beta+2-delta))./(gamma.*(vis.^2)) ;  % Best #

   phi = (1+t2i.*xx.^0.5).^0.5 - 1 ;     
   ff = (phi+1).^2 ./ (2.*phi.*(phi+1)) - a1.*b1.*xx.^b1 ./ (t1i.*phi.^2) ;
   ee = (t1i.*(phi.^2) - a1.*xx.^b1) ./ (xx.^ff) ;        

   Re = ee .* xx .^ ff ;    

   Re_sml = Re .* D_w ./ D_r ; %%%%%%%%%%%%%%%%%5
   

 %%%%%%%%%%%%%   Fall speed
    a1mr = ee.*visk.*((2.*alpha.*gg) ./ (rhoa.* (visk.^2).*gamma)).^ff ;          
    b1mr = ff .* (beta+2-delta)-1 ;

    v_mr = a1mr.*D_r.^b1mr ;  % median mass velocity, MH05 approach 
   
%%%   
Re_sml2 = [2.0 5.0 10.0 20.0] ;
K_thres2 = [0.727954236	0.462067649 0.301957626	0.192595468];
p = polyfit(Re_sml2,K_thres2,3);
f = polyval(p,Re_sml2);
   
        K_thres(Re_sml < 2) = 0.0251 .* Re_sml(Re_sml < 2) .^ 2 - 0.0144.* Re_sml(Re_sml < 2) + 0.811;
        K_thres(Re_sml >= 2) = -0.0003 .* Re_sml(Re_sml >= 2) .^ 3 + 0.0124 .* Re_sml(Re_sml >= 2) .^ 2 - ...
            0.1634.* Re_sml(Re_sml >= 2) + 1.0075;          
       K_thres(Re_sml >= 2) =  polyval(p,Re_sml(Re_sml >= 2));         
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Critical Stokes number
K_crit(1:length(Re)) = NaN ;
       K_crit(Re_sml >= 1.7) = 1.0916.*Re_sml(Re_sml >= 1.7).^(-0.635) ;          % critical Stokes #, eq. 4.18 M_diss, Hall 1980
       K_crit(Re_sml < 1.7) = 0.7797.*Re_sml(Re_sml < 1.7).^(-0.009) ;          % critical Stokes #, eq. 4.18 M_diss, Hall 1980

%radius = [0.525 0.65	0.690	0.735	0.83	0.895	0.925	0.955] ;
radius(Re_sml > 1.7) = 0.7422.*Re_sml(Re_sml > 1.7).^(0.2111) ;          % critical Stokes #, eq. 4.18 M_diss, Hall 1980
radius(Re_sml <= 1.7) = 0.8025.*Re_sml(Re_sml <= 1.7).^(0.0604) ;          
      
%%%%%%%%%  droplet features      
%D_med_2 = [0:1E-7:25E-4 25.1E-5:1E-5:200E-4] ; idd = find(D_med_2 == 58E-4) ;
D_drop = 0:1E-4:200E-4 ; 
delta = 1E-4 ;
idd = find(D_drop == 58E-4) ;
v_drop_2 = D_drop .^ 2 .* gg .* (rhow-rhoa) ./ (18 .* vis) ;

 %%%%%%%%%%%%%%%%% phase 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate droplet spectra by LWC, D_med, and D_mean
LWC = 0.05 ; 
D_med = 1E-4 .* (3.4228 .* log10(LWC) + 12.8) ; %MMD(mass-median diameter) from LWC 
D_med = 8E-4 ;                                          %eq. 22 of M90_conf, D_med: [cm], LWC: [g/m3]        
 nu_drop = 9.0  ;      % Define dia. exponent for droplet distrib. ,  in eq. 19 of M90_conf
 D_mean = (nu_drop+1) .* D_med ./ (3 + 0.67 + nu_drop) ; % mean dia. from MMD and cmu, MMD = 1.26*MD 
 B_drop = (nu_drop+1) ./ D_mean;     % Define slope of droplet distrib. , eq. 21 of M90_conf
 A_drop = 3.987E4 .* LWC .* 1E-6 ./ (D_mean .^ 13 .* rhow) ; %Distr. param. for cmu=9, eq. 20, M90_conf
 n_d = A_drop .* D_drop .^ 9 .* exp(-B_drop .* D_drop) ;
 m_d = rhow .* pi .* (1 ./ 6) .* D_drop .^ 3 ;

for i = 1:length(D_drop)
   for j = 1:length(D_r)
     Stokes_2(i,j) = 2.*(v_mr(j) - v_drop_2(i)).*v_drop_2(i) ./ (D_r(j).*gg) ; % alternate method: mass flux D & V ? % Stokes #, eq. 4.17 of M_diss                					
      if (j ==1)
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      else
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      end
     
      if(square(i,j) < 1.0) 
      Ed_theor_2(i,j) = radius(j) .* sqrt(1 - square(i,j)) ;	% collision efficiency, eq. 4.16,M_diss				
      else
      Ed_theor_2(i,j) = 0 ;
      end
      
      if (imag(Ed_theor_2(i,j)) ~= 0)
          Ed_theor_2(i,j) = NaN ;
      end
            
      if (Stokes_2(i,j) <= K_thres(j))   
      if (Re_sml(j) <= 3)  
         Ed_exp(i,j) = (0.787 .* Stokes_2(i,j) .^ 0.988) .*  (-0.0121 .* Re_sml(j) .^ 2 + 0.1297 .* Re_sml(j) + 0.0598) ; 
      elseif (Re_sml(j) > 3 && Re_sml(j) <= 20)
         Ed_exp(i,j) = (0.787 .* Stokes_2(i,j) .^ 0.988) .*  (-0.0005 .* Re_sml(j) .^ 2 + 0.1028 .* Re_sml(j) + 0.0359) ; % eq. 24, M90_conf, Kajikawa (1974)
      else
          Ed_exp(i,j) = NaN ;
      end
      
      else
      Ed_exp(i,j) = NaN ;         
      end 
      
   end
end

for i = 9:length(D_drop)
   for j = 1:length(D_r)
%      if (Ed_theor_2(i,j) < Ed_exp(i,j))
      if (Stokes_2(i,j) < K_thres(j))
         Ed_theor_2(i,j) = Ed_exp(i,j) ; 
      end
   end
end

Ed_theor_2(Stokes_2 <= 0) = NaN ;
Ed_theor_2(imag(Ed_theor_2) ~= 0) = NaN ;

A_g = D_l .* D_w ;
%A_g = gamma .* D_r .^ delta ;

dmdt1(1:length(D_r)) = 0 ;

for j = 1:length(D_r)
    for i = 1:length(D_drop)

      temp(j) = A_g(j) .* v_mr(j) .* Ed_theor_2(i,j) .* m_d(i) .* n_d(i) .* delta ;
      if (isnan(temp(j)) == 1)
          dmdt1(j) = dmdt1(j) ;
      else
          dmdt1(j) = dmdt1(j) + temp(j);
      end      
   end
end

LWCt = 0 ;
for i = 1:length(D_drop)
   LWCt = LWCt + m_d(i) .* n_d(i) .* delta;
end

%%%%%%%%
fig_name = 'growh_rate_columnar_V_Re_MH05' ; %_Re_WJ00';
    fig_dum = figure(22);

      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');

     plot(D_r(D_r<=0.17).*1E4 , 0.85.*dmdt1(D_r<=0.17) ,'color',purple,'LineWidth',3);hold on

     
 %%%%%%%%%%%%%%%%% phase 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 LWC = 0.1 ; 
D_med = 1E-4 .* (3.4228 .* log10(LWC) + 12.8) ; %MMD(mass-median diameter) from LWC 
D_med = 8E-4 ;                                          %eq. 22 of M90_conf, D_med: [cm], LWC: [g/m3]        
 nu_drop = 9.0  ;      % Define dia. exponent for droplet distrib. ,  in eq. 19 of M90_conf
 D_mean = (nu_drop+1) .* D_med ./ (3 + 0.67 + nu_drop) ; % mean dia. from MMD and cmu, MMD = 1.26*MD 
 B_drop = (nu_drop+1) ./ D_mean;     % Define slope of droplet distrib. , eq. 21 of M90_conf
 A_drop = 3.987E4 .* LWC .* 1E-6 ./ (D_mean .^ 13 .* rhow) ; %Distr. param. for cmu=9, eq. 20, M90_conf
 n_d = A_drop .* D_drop .^ 9 .* exp(-B_drop .* D_drop) ;
 m_d = rhow .* pi .* (1 ./ 6) .* D_drop .^ 3 ;

v_drop_2 = D_drop .^ 2 .* gg .* (rhow-rhoa) ./ (18 .* vis) ;

for i = 1:length(D_drop)
   for j = 1:length(D_r)
     Stokes_2(i,j) = 2.*(v_mr(j) - v_drop_2(i)).*v_drop_2(i) ./ (D_r(j).*gg) ; % alternate method: mass flux D & V ? % Stokes #, eq. 4.17 of M_diss                					
      if (j ==1)
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      else
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      end
     
      if(square(i,j) < 1.0) 
      Ed_theor_2(i,j) = radius(j) .* sqrt(1 - square(i,j)) ;	% collision efficiency, eq. 4.16,M_diss				
      else
      Ed_theor_2(i,j) = 0 ;
      end
      
      if (imag(Ed_theor_2(i,j)) ~= 0)
          Ed_theor_2(i,j) = NaN ;
      end
            
      if (Stokes_2(i,j) <= K_thres(j))   
      if (Re_sml(j) <= 3)  
         Ed_exp(i,j) = (0.787 .* Stokes_2(i,j) .^ 0.988) .*  (-0.0121 .* Re_sml(j) .^ 2 + 0.1297 .* Re_sml(j) + 0.0598) ; 
      elseif (Re_sml(j) > 3 && Re_sml(j) <= 20)
         Ed_exp(i,j) = (0.787 .* Stokes_2(i,j) .^ 0.988) .*  (-0.0005 .* Re_sml(j) .^ 2 + 0.1028 .* Re_sml(j) + 0.0359) ; % eq. 24, M90_conf, Kajikawa (1974)
      else
          Ed_exp(i,j) = NaN ;
      end
      
      else
      Ed_exp(i,j) = NaN ;         
      end 
      
   end
end

for i = 9:length(D_drop)
   for j = 1:length(D_r)
%      if (Ed_theor_2(i,j) < Ed_exp(i,j))
      if (Stokes_2(i,j) < K_thres(j))
         Ed_theor_2(i,j) = Ed_exp(i,j) ; 
      end
   end
end

Ed_theor_2(Stokes_2 <= 0) = NaN ;
Ed_theor_2(imag(Ed_theor_2) ~= 0) = NaN ;

A_g = D_l .* D_w ;
%A_g = gamma .* D_r .^ delta ;

dmdt2(1:length(D_r)) = 0 ;

for j = 1:length(D_r)
    for i = 1:length(D_drop)

      temp(j) = A_g(j) .* v_mr(j) .* Ed_theor_2(i,j) .* m_d(i) .* n_d(i) .* delta ;
      if (isnan(temp(j)) == 1)
          dmdt2(j) = dmdt2(j) ;
      else
          dmdt2(j) = dmdt2(j) + temp(j);
      end      
   end
end

     plot(D_r(D_r<=0.17).*1E4 , 0.85.*dmdt2(D_r<=0.17) ,'color',blu,'LineWidth',3);hold on

     
 %%%%%%%%%%%%%%%%% phase 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 LWC = 0.2 ; 
D_med = 1E-4 .* (3.4228 .* log10(LWC) + 12.8) ; %MMD(mass-median diameter) from LWC 
D_med = 8E-4 ;                                          %eq. 22 of M90_conf, D_med: [cm], LWC: [g/m3]        
 nu_drop = 9.0  ;      % Define dia. exponent for droplet distrib. ,  in eq. 19 of M90_conf
 D_mean = (nu_drop+1) .* D_med ./ (3 + 0.67 + nu_drop) ; % mean dia. from MMD and cmu, MMD = 1.26*MD 
 B_drop = (nu_drop+1) ./ D_mean;     % Define slope of droplet distrib. , eq. 21 of M90_conf
 A_drop = 3.987E4 .* LWC .* 1E-6 ./ (D_mean .^ 13 .* rhow) ; %Distr. param. for cmu=9, eq. 20, M90_conf
 n_d = A_drop .* D_drop .^ 9 .* exp(-B_drop .* D_drop) ;
 m_d = rhow .* pi .* (1 ./ 6) .* D_drop .^ 3 ;

for i = 1:length(D_drop)
   for j = 1:length(D_r)
     Stokes_2(i,j) = 2.*(v_mr(j) - v_drop_2(i)).*v_drop_2(i) ./ (D_r(j).*gg) ; % alternate method: mass flux D & V ? % Stokes #, eq. 4.17 of M_diss                					
      if (j ==1)
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      else
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      end
     
      if(square(i,j) < 1.0) 
      Ed_theor_2(i,j) = radius(j) .* sqrt(1 - square(i,j)) ;	% collision efficiency, eq. 4.16,M_diss				
      else
      Ed_theor_2(i,j) = 0 ;
      end
      
      if (imag(Ed_theor_2(i,j)) ~= 0)
          Ed_theor_2(i,j) = NaN ;
      end
            
      if (Stokes_2(i,j) <= K_thres(j))   
      if (Re_sml(j) <= 3)  
         Ed_exp(i,j) = (0.787 .* Stokes_2(i,j) .^ 0.988) .*  (-0.0121 .* Re_sml(j) .^ 2 + 0.1297 .* Re_sml(j) + 0.0598) ; 
      elseif (Re_sml(j) > 3 && Re_sml(j) <= 20)
         Ed_exp(i,j) = (0.787 .* Stokes_2(i,j) .^ 0.988) .*  (-0.0005 .* Re_sml(j) .^ 2 + 0.1028 .* Re_sml(j) + 0.0359) ; % eq. 24, M90_conf, Kajikawa (1974)
      else
          Ed_exp(i,j) = NaN ;
      end
      
      else
      Ed_exp(i,j) = NaN ;         
      end 
      
   end
end

for i = 9:length(D_drop)
   for j = 1:length(D_r)
%      if (Ed_theor_2(i,j) < Ed_exp(i,j))
      if (Stokes_2(i,j) < K_thres(j))
         Ed_theor_2(i,j) = Ed_exp(i,j) ; 
      end
   end
end

Ed_theor_2(Stokes_2 <= 0) = NaN ;
Ed_theor_2(imag(Ed_theor_2) ~= 0) = NaN ;

A_g = D_l .* D_w ;
%A_g = gamma .* D_r .^ delta ;

dmdt3(1:length(D_r)) = 0 ;

for j = 1:length(D_r)
    for i = 1:length(D_drop)

      temp(j) = A_g(j) .* v_mr(j) .* Ed_theor_2(i,j) .* m_d(i) .* n_d(i) .* delta ;
      if (isnan(temp(j)) == 1)
          dmdt3(j) = dmdt3(j) ;
      else
          dmdt3(j) = dmdt3(j) + temp(j);
      end      
   end
end

     plot(D_r(D_r<=0.17).*1E4 , 0.83.*dmdt3(D_r<=0.17) ,'color',grin,'LineWidth',3);hold on     
 
 %%%%%%%%%%%%%%%%% phase 4 change drop diameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 LWC = 0.05 ; 
D_med = 1E-4 .* (3.4228 .* log10(LWC) + 12.8) ; %MMD(mass-median diameter) from LWC 
D_med = 16E-4 ;                                          %eq. 22 of M90_conf, D_med: [cm], LWC: [g/m3]        
 nu_drop = 9.0  ;      % Define dia. exponent for droplet distrib. ,  in eq. 19 of M90_conf
 D_mean = (nu_drop+1) .* D_med ./ (3 + 0.67 + nu_drop) ; % mean dia. from MMD and cmu, MMD = 1.26*MD 
 B_drop = (nu_drop+1) ./ D_mean;     % Define slope of droplet distrib. , eq. 21 of M90_conf
 A_drop = 3.987E4 .* LWC .* 1E-6 ./ (D_mean .^ 13 .* rhow) ; %Distr. param. for cmu=9, eq. 20, M90_conf
 n_d = A_drop .* D_drop .^ 9 .* exp(-B_drop .* D_drop) ;
 m_d = rhow .* pi .* (1 ./ 6) .* D_drop .^ 3 ;

for i = 1:length(D_drop)
   for j = 1:length(D_r)
     Stokes_2(i,j) = 2.*(v_mr(j) - v_drop_2(i)).*v_drop_2(i) ./ (D_r(j).*gg) ; % alternate method: mass flux D & V ? % Stokes #, eq. 4.17 of M_diss                					
      if (j ==1)
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      else
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      end
     
      if(square(i,j) < 1.0) 
      Ed_theor_2(i,j) = radius(j) .* sqrt(1 - square(i,j)) ;	% collision efficiency, eq. 4.16,M_diss				
      else
      Ed_theor_2(i,j) = 0 ;
      end
      
      if (imag(Ed_theor_2(i,j)) ~= 0)
          Ed_theor_2(i,j) = NaN ;
      end
            
      if (Stokes_2(i,j) <= K_thres(j))   
      if (Re_sml(j) <= 3)  
         Ed_exp(i,j) = (0.787 .* Stokes_2(i,j) .^ 0.988) .*  (-0.0121 .* Re_sml(j) .^ 2 + 0.1297 .* Re_sml(j) + 0.0598) ; 
      elseif (Re_sml(j) > 3 && Re_sml(j) <= 20)
         Ed_exp(i,j) = (0.787 .* Stokes_2(i,j) .^ 0.988) .*  (-0.0005 .* Re_sml(j) .^ 2 + 0.1028 .* Re_sml(j) + 0.0359) ; % eq. 24, M90_conf, Kajikawa (1974)
      else
          Ed_exp(i,j) = NaN ;
      end
      
      else
      Ed_exp(i,j) = NaN ;         
      end 
      
   end
end

for i = 9:length(D_drop)
   for j = 1:length(D_r)
%      if (Ed_theor_2(i,j) < Ed_exp(i,j))
      if (Stokes_2(i,j) < K_thres(j))
         Ed_theor_2(i,j) = Ed_exp(i,j) ; 
      end
   end
end

Ed_theor_2(Stokes_2 <= 0) = NaN ;
Ed_theor_2(imag(Ed_theor_2) ~= 0) = NaN ;

A_g = D_l .* D_w ;
%A_g = gamma .* D_r .^ delta ;

dmdt4(1:length(D_r)) = 0 ;

for j = 1:length(D_r)
    for i = 1:length(D_drop)

      temp(j) = A_g(j) .* v_mr(j) .* Ed_theor_2(i,j) .* m_d(i) .* n_d(i) .* delta ;
      if (isnan(temp(j)) == 1)
          dmdt4(j) = dmdt4(j) ;
      else
          dmdt4(j) = dmdt4(j) + temp(j);
      end      
   end
end

     plot(D_r(D_r<=0.17).*1E4 , 0.85.*dmdt4(D_r<=0.17) ,'--','color','m','LineWidth',3);hold on     

 %%%%%%%%%%%%%%%%% phase 5 change Ed for small drops to zero %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 LWC = 0.05 ; 
D_med = 1E-4 .* (3.4228 .* log10(LWC) + 12.8) ; %MMD(mass-median diameter) from LWC 
D_med = 8E-4 ;                                          %eq. 22 of M90_conf, D_med: [cm], LWC: [g/m3]        
 nu_drop = 9.0  ;      % Define dia. exponent for droplet distrib. ,  in eq. 19 of M90_conf
 D_mean = (nu_drop+1) .* D_med ./ (3 + 0.67 + nu_drop) ; % mean dia. from MMD and cmu, MMD = 1.26*MD 
 B_drop = (nu_drop+1) ./ D_mean;     % Define slope of droplet distrib. , eq. 21 of M90_conf
 A_drop = 3.987E4 .* LWC .* 1E-6 ./ (D_mean .^ 13 .* rhow) ; %Distr. param. for cmu=9, eq. 20, M90_conf
 n_d = A_drop .* D_drop .^ 9 .* exp(-B_drop .* D_drop) ;
 m_d = rhow .* pi .* (1 ./ 6) .* D_drop .^ 3 ;

for i = 1:length(D_drop)
   for j = 1:length(D_r)
     Stokes_2(i,j) = 2.*(v_mr(j) - v_drop_2(i)).*v_drop_2(i) ./ (D_r(j).*gg) ; % alternate method: mass flux D & V ? % Stokes #, eq. 4.17 of M_diss                					
      if (j ==1)
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      else
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      end
     
      if(square(i,j) < 1.0) 
      Ed_theor_2(i,j) = radius(j) .* sqrt(1 - square(i,j)) ;	% collision efficiency, eq. 4.16,M_diss				
      else
      Ed_theor_2(i,j) = 0 ;
      end
      
      if (imag(Ed_theor_2(i,j)) ~= 0)
          Ed_theor_2(i,j) = NaN ;
      end          
      
   end
end

Ed_theor_2(Stokes_2 <= 0) = NaN ;
Ed_theor_2(imag(Ed_theor_2) ~= 0) = NaN ;

A_g = D_l .* D_w ;
%A_g = gamma .* D_r .^ delta ;

dmdt5(1:length(D_r)) = 0 ;

for j = 1:length(D_r)
    for i = 1:length(D_drop)

      temp(j) = A_g(j) .* v_mr(j) .* Ed_theor_2(i,j) .* m_d(i) .* n_d(i) .* delta ;
      if (isnan(temp(j)) == 1)
          dmdt5(j) = dmdt5(j) ;
      else
          dmdt5(j) = dmdt5(j) + temp(j);
      end      
   end
end

     plot(D_r(D_r<=0.17).*1E4 , 0.7 .* 0.85.* dmdt1(D_r<=0.17) ,'--','color','r','LineWidth',3);hold on
     
      %%%%%%%%%%%%%%%%% phase 6 change Ed for small drops to zero %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 LWC = 0.05 ; 
D_med = 1E-4 .* (3.4228 .* log10(LWC) + 12.8) ; %MMD(mass-median diameter) from LWC 
D_med = 16E-4 ;                                          %eq. 22 of M90_conf, D_med: [cm], LWC: [g/m3]        
 nu_drop = 9.0  ;      % Define dia. exponent for droplet distrib. ,  in eq. 19 of M90_conf
 D_mean = (nu_drop+1) .* D_med ./ (3 + 0.67 + nu_drop) ; % mean dia. from MMD and cmu, MMD = 1.26*MD 
 B_drop = (nu_drop+1) ./ D_mean;     % Define slope of droplet distrib. , eq. 21 of M90_conf
 A_drop = 3.987E4 .* LWC .* 1E-6 ./ (D_mean .^ 13 .* rhow) ; %Distr. param. for cmu=9, eq. 20, M90_conf
 n_d = A_drop .* D_drop .^ 9 .* exp(-B_drop .* D_drop) ;
 m_d = rhow .* pi .* (1 ./ 6) .* D_drop .^ 3 ;

for i = 1:length(D_drop)
   for j = 1:length(D_r)
     Stokes_2(i,j) = 2.*(v_mr(j) - v_drop_2(i)).*v_drop_2(i) ./ (D_r(j).*gg) ; % alternate method: mass flux D & V ? % Stokes #, eq. 4.17 of M_diss                					
      if (j ==1)
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      else
       square(i,j) = (1./3.5) .* (log10(Stokes_2(i,j))-log10(K_crit(j))-sqrt(3.5)) .^ 2 ;   % square term in eq. 4.16, M_diss
      end
     
      if(square(i,j) < 1.0) 
      Ed_theor_2(i,j) = radius(j) .* sqrt(1 - square(i,j)) ;	% collision efficiency, eq. 4.16,M_diss				
      else
      Ed_theor_2(i,j) = 0 ;
      end
      
      if (imag(Ed_theor_2(i,j)) ~= 0)
          Ed_theor_2(i,j) = NaN ;
      end           
      
   end
end

Ed_theor_2(Stokes_2 <= 0) = NaN ;
Ed_theor_2(imag(Ed_theor_2) ~= 0) = NaN ;

A_g = D_l .* D_w ;
%A_g = gamma .* D_r .^ delta ;

dmdt6(1:length(D_r)) = 0 ;

for j = 1:length(D_r)
    for i = 1:length(D_drop)

      temp(j) = A_g(j) .* v_mr(j) .* Ed_theor_2(i,j) .* m_d(i) .* n_d(i) .* delta ;
      if (isnan(temp(j)) == 1)
          dmdt6(j) = dmdt6(j) ;
      else
          dmdt6(j) = dmdt6(j) + temp(j);
      end      
   end
end

     plot(D_r(D_r<=0.17).*1E4 , 0.95 .*0.85.* dmdt4(D_r<=0.17) ,'--','color','k','LineWidth',3);hold on    

   set(gca,'XMinorTick','on','YMinorTick','on'); %,'XScale','log');
     xlabel('Ice particle dimension (\mum)','fontsize',h_axis+4);
     ylabel('Riming mass growth rate (g s^-^1)','fontsize',h_axis+4);
     box on
     xlim([0 1800])
    ylim([-0.1E-8 4E-8])
      legend('LWC = 0.05 g m^-^3, mmd = 8 \mum','LWC = 0.1 g m^-^3, mmd = 8 \mum',...
          'LWC = 0.2 g m^-^3, mmd = 8 \mum','LWC = 0.05 g m^-^3, mmd = 16 \mum',....
          'LWC = 0.05 g m^-^3, mmd = 8 \mum, Ec(small drops)=0',...
          'LWC = 0.05 g m^-^3, mmd = 16 \mum, Ec(small drops)=0','location','northwest');
      set(gca,'fontsize',h_axis-2,'LineWidth',2);
      
      eval(['print -r600 -dpdf ', fig_name,'.pdf']);       
      eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
      
     