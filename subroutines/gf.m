% the green function of spherically layered media
% use Octave's library to calculate associated Legendre functions
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

function field_value = gf(cal,u,e,a,r,theta,phi,r_p,theta_p,phi_p,k0,precision)
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
field_value = 0;
delta = 0;
delta_old = 1e10;

switch cal
case {1} % calculate H_r

  do	
    sp = legendre(n,x,'sch');
    sp_p = legendre(n,x_p,'sch');
    % the definition of sp^0_0 in octave is different with mine.
    sp(1) = sqrt(2)*sp(1); 
    sp_p(1) = sqrt(2)*sp_p(1);
    % calculate delta which is a summation over m
    % m=0 case
    Dsp_p = sqrt( n*(n+1) ) * -sp_p(2) ;
    delta = sp(1) * Dsp_p ;
    % m=1:n case
    for m=1:n
      Dsp_p = sqrt( (n+m)*(n-m+1) )*sp_p(m) - m*cot(theta_p)*sp_p(m+1);
      delta = delta + 2 * sp(m+1) * Dsp_p * cos(m*dphi) ;
    end
    [F dF_rp dF_r d2F] = F_tetm('te',u,e,a,r,r_p,k0,n);
    delta = (n+0.5)*F*delta;
    field_value = field_value + delta; % add to total field
    
    if ( abs( delta/field_value )<precision && abs( delta_old/field_value )<precision ) % reach the prcision
      printf('H_r: converged after %d iterations. \n', n);
      break
    else % otherwise, update for n+1 
      n = n+1; 
      delta_old = delta; 
    end
  until (0)
  field_value = field_value * u(jj)/u(ii)*kj/(4i*pi*r) ; 

case {2} % calculate E_r

  do	
    sp = legendre(n,x,'sch');
    sp_p = legendre(n,x_p,'sch');
    % calculate delta which is a summation over m
    for m=1:n
      delta = delta + m*sp(m+1)*sp_p(m+1)*sin(m*dphi) ;
    end
    [F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_p,k0,n);
    delta = (2*n+1)*(F/r_p+dF_rp)*delta;
    field_value = field_value + delta; % add to total field
    
    if ( abs( delta/field_value )<precision && abs( delta_old/field_value )<precision ) % reach the prcision
      printf('E_r: converged after %d iterations. \n', n);
      break
    else % otherwise, update for n+1 
      n = n+1; 
      delta_old = delta; 
      delta = 0;
    end
  until (0)
  field_value = field_value * -kj/(4*pi*sin(theta_p)) /(w*e(ii)*e0*r); 

end



