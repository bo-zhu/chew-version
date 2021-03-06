% the green function of a PHI direction dipole in spherically layered media
% use recurrence formulas to calculate associated Legendre functions
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

function field_value = gf_phi_SPrecursion(cal,u,e,a,r,theta,phi,r_p,theta_p,phi_p,k0,precision)
% cal: 1 H_r; 2 E_r
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
sp_p_init = sin(theta_p); % initial value  
sp = sp_init;
sp_p = sp_p_init;
sp_old = 0; % sp^2_1
sp_p_old = 0;  

field_value = 0;
delta = 0;
delta_old = 1e10;

switch cal
case {1} % calculate H_r

  do	
    for m=n:-1:0
      % calculate delta which is a summation over m
      sp_p_new = legendre_decrease_m('sch',x_p,sp_p,sp_p_old,n,m); % sp_p^m-1_n
      Dsp_p = sqrt( (n+m)*(n-m+1) ) * sp_p_new - m*cot(theta_p)*sp_p;
      delta = delta + (1+min(1,m)) * sp * Dsp_p * cos(m*dphi) ;
      % update sp and sp_p for m-1
      sp_new = legendre_decrease_m('sch',x,sp,sp_old,n,m); 
      sp_old = sp;
      sp = sp_new;
      sp_p_old = sp_p;
      sp_p = sp_p_new;
    end
    [F dF_rp dF_r d2F] = F_tetm('te',u,e,a,r,r_p,k0,n);
    delta = (n+0.5)*F*delta;
    field_value = field_value + delta; % add to total field
    
    if ( abs( delta/field_value )<precision && abs( delta_old/field_value )<precision ) % reach the precision
      printf('H_r: converged after %d iterations. \n', n);
      break
    else % otherwise, update for n+1 
      n = n+1; 
      b1 = b1+2;
      b2 = b2+2;
      sp_init = sp_init*sqrt(b1/b2)*sin(theta);
      sp_p_init = sp_p_init*sqrt(b1/b2)*sin(theta_p);
      sp = sp_init;
      sp_p = sp_p_init;
      sp_old = 0;
      sp_p_old = 0;
      delta_old = delta; 
      delta = 0;
    end
  until (0)
  field_value = field_value * u(jj)/u(ii)*kj/(4i*pi*r) ; 

case {2} % calculate E_r

  do	
    for m=n:-1:1
      % calculate delta which is a summation over m
      delta = delta + m*sp*sp_p*sin(m*dphi) ;
      % update sp and sp_p for m-1
      sp_new = legendre_decrease_m('sch',x,sp,sp_old,n,m); 
      sp_old = sp;
      sp = sp_new;
      sp_p_new = legendre_decrease_m('sch',x_p,sp_p,sp_p_old,n,m); 
      sp_p_old = sp_p;
      sp_p = sp_p_new;
    end
    [F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_p,k0,n);
    delta = (2*n+1)*(F/r_p+dF_rp)*delta;
    field_value = field_value + delta; % add to total field
    
    if ( abs( delta/field_value )<precision && abs( delta_old/field_value )<precision ) % reach the precision
      printf('E_r: converged after %d iterations. \n', n);
      break
    else % otherwise, update for n+1 
      n = n+1; 
      b1 = b1+2;
      b2 = b2+2;
      sp_init = sp_init*sqrt(b1/b2)*sin(theta);
      sp_p_init = sp_p_init*sqrt(b1/b2)*sin(theta_p);
      sp = sp_init;
      sp_p = sp_p_init;
      sp_old = 0;
      sp_p_old = 0;
      delta_old = delta; 
      delta = 0;
    end
  until (0)
  field_value = field_value * -kj/(4*pi*sin(theta_p)) /(w*e(ii)*e0*r); 

end


