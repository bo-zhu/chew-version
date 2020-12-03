% the green function of an R direction dipole in spherically layered media
% use recurrence formulas to calculate associated Legendre functions
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

function field_value = gf_r_SPrecursion(cal,u,e,a,r,theta,phi,r_p,theta_p,phi_p,k0,precision)
% cal: 2 E_r
% u,e : the relative medium parameters of each shell.
% a : the radius of each interface of neighboring shells.
% r,theta,phi : the observation point's coordinate. 
% r_p,theta_p,phi_p : the source point's coordinate. 
% k0 : vacuum wave number.
% precision : the precision to stop the iteration.


%%%%%%%% physical constants %%%%%%%%%%
e0 = 8.854e-12;
u0 = 4*pi*1e-7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = find(a>r);
ii = idx(1);
idx = find(a>r_p);
jj = idx(1);

w = k0/sqrt(e0*u0);
f = w/2/pi;
kj = k0*sqrt(e(jj)*u(jj));
ki = k0*sqrt(e(ii)*u(ii));
	
x = cos(theta);
x_p = cos(theta_p);
dphi = phi-phi_p;

n = 1;
b1 = 1;
b2 = 2;
sp_init = sin(theta); % initial value: sp^1_1
sp = sp_init;
sp_old = 0; % sp^2_1
sp_p_init = sin(theta_p); % initial value: sp_p^1_1  
sp_p = sp_p_init;
sp_p_old = 0; % sp_p^2_1
sp_p_ = legendre_decrease_m('sch',x_p,sp_p,sp_p_old,n,n); % sp_p^0_1

field_value = 0;
delta = 0;
delta_old = 1e10;

switch cal
case {2} % calculate E_r

  sin_t = sin(theta);
  sin_p = sin(theta_p);
  cot_p = cot(theta_p);
  sin2_p = sin_p^2;
  csc2_p = csc(theta_p)^2;
  do	
    % calculate delta which is a summation over m
    %%%%%% m=n:-1:0 case %%%%%%
    for m=n:-1:0
      if n==1&&m==0
	sp_p__ = 0;
      else
        sp_p__ = legendre_decrease_m('sch',x_p,sp_p_,sp_p,n,m-1); % sp_p^m-2_n
      end
      Dsp_p = sqrt( (n+m)*(n-m+1) )*sp_p_ - m*cot_p*sp_p;
      Dsp_p_ = sqrt( (n+m-1)*(n-m+2) )*sp_p__ - (m-1)*cot_p*sp_p_;
      DDsp_p = sqrt( (n+m)*(n-m+1) )*Dsp_p_ + m*csc2_p*sp_p - m*cot_p*Dsp_p;

      delta = delta + (1+min(1,m))*sp * ( cot_p*Dsp_p + DDsp_p - m^2*sp_p/sin2_p ) * cos(m*dphi) ;
      % update sp_p 
      sp_p = sp_p_;
      sp_p_ = sp_p__;
      % update sp 
      sp_ = legendre_decrease_m('sch',x,sp,sp_old,n,m); % sp^m-1_n
      sp_old = sp;
      sp = sp_;
    end

    [F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_p,k0,n);
    delta = (n+0.5)*F*delta;
    field_value = field_value + delta; % add to total field
    if ( abs( delta/field_value )<precision && abs( delta_old/field_value )<precision ) % reach the precision
      printf('E_r: converged after %d iterations. \n', n);
      break
    else % otherwise, update for n+1 
      n = n+1; 
      b1 = b1+2;
      b2 = b2+2;
      sp_init = sp_init*sqrt(b1/b2)*sin_t;
      sp = sp_init;
      sp_old = 0;
      sp_p_init = sp_p_init*sqrt(b1/b2)*sin_p;
      sp_p = sp_p_init;
      sp_p_old = 0;
      sp_p_ = legendre_decrease_m('sch',x_p,sp_p,sp_p_old,n,n); % sp_p^n-1_n
      delta_old = delta; 
      delta = 0;
    end
  until (0)

  field_value =  kj/(4*pi*w*e(ii)*e0*r*r_p) * field_value ; 


end



