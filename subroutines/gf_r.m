% the green function of an R direction dipole in spherically layered media
% use Octave's library to calculate associated Legendre functions
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

function field_value = gf_r(cal,u,e,a,r,theta,phi,r_p,theta_p,phi_p,k0,precision)
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

switch cal
case {2} % calculate E_r

  SUM = 0;
  n = 1;
  delta = 0;
  delta_old = 1e10;
  cot_p = cot(theta_p);
  sin2_p = sin(theta_p)^2;
  csc2_p = csc(theta_p)^2;
  do	
    sp = legendre(n,x,'sch');
    sp_p = legendre(n,x_p,'sch');
    % the definition of sp^0_0 in octave is different with mine.
    sp(1) = sqrt(2)*sp(1); 
    sp_p(1) = sqrt(2)*sp_p(1);

    % calculate delta which is a summation over m
    %%%%%% m=0 case %%%%%%
    Dsp_p = sqrt( n*(n+1) ) * -sp_p(2) ;
    if n==1
      DDsp_p = sqrt( n*(n+1) ) * ( -cot_p*sp_p(2) );
    else
      DDsp_p = sqrt( n*(n+1) ) * ( sqrt( (n-1)*(n+2) )*sp_p(3) - cot_p*sp_p(2) );
    end
    delta = sp(1) * ( cot_p*Dsp_p + DDsp_p );
    %%%%%% m=1 case %%%%%%
    Dsp_p = sqrt( (n+1)*n )*sp_p(1) - cot_p*sp_p(2);
    DDsp_p = ( 2*csc2_p-1-n*(n+1) )*sp_p(2) - cot_p*sp_p(1)*sqrt(n*(n+1));
    delta = delta + 2*sp(2) * ( cot_p*Dsp_p+DDsp_p-sp_p(2)/sin2_p ) * cos(dphi) ;
    %%%%%% m=2:n case %%%%%%
    for m=2:n
      Dsp_p = sqrt( (n+m)*(n-m+1) )*sp_p(m) - m*cot_p*sp_p(m+1);
      Dsp_p_ = sqrt( (n+m-1)*(n-m+2) )*sp_p(m-1) - (m-1)*cot_p*sp_p(m);
      DDsp_p = sqrt( (n+m)*(n-m+1) )*Dsp_p_ + m*csc2_p*sp_p(m+1) - m*cot_p*Dsp_p;
      delta = delta + 2*sp(m+1) * ( cot_p*Dsp_p + DDsp_p - m^2*sp_p(m+1)/sin2_p ) * cos(m*dphi) ;
    end

    [F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_p,k0,n);
    delta = (n+0.5)*F*delta;
    SUM = SUM + delta; % add to total field
    if ( abs( delta/SUM )<precision && abs( delta_old/SUM )<precision ) % reach the precision
      printf('E_r: converged after %d iterations. \n', n);
      break
    else % otherwise, update for n+1 
      n = n+1; 
      delta_old = delta; 
    end
  until (0)

  field_value =  kj/(4*pi*w*e(ii)*e0*r*r_p) * SUM; 


end



