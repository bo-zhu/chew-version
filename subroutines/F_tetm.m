% calculate Fn(r,r'), dFn/dr', dFn/dr, and d^2Fn/dr'dr.
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

function [F dF_rp dF_r d2F] = F_tetm(mode,u,e,a,r,rp,k0,n)
% u,e : the relative medium parameters of each shell.
% a : the radius of each interface of neighboring shells.
% r : the observation point's radius. 
% rp : the radius of the charge's orbit.
% n : the degree of the spherical Bessel function. 


idx = find(a>r);
ii = idx(1);
u_out_ii = u(ii:end);
e_out_ii = e(ii:end);
a_out_ii = a(ii:end);
u_std_ii = u(1:ii); 
e_std_ii = e(1:ii);
a_std_ii = a(1:ii-1);

idx = find(a>rp);
jj = idx(1);
u_out_jj = u(jj:end);
e_out_jj = e(jj:end);
a_out_jj = a(jj:end);
u_std_jj = u(1:jj); 
e_std_jj = e(1:jj);
a_std_jj = a(1:jj-1);

kj = k0*sqrt( u(jj)*e(jj) );

%%%%%%%% calculate parameters relevant to aj, aj__, ai and ai__ %%%%%%%%%
kj_aj__ = kj*a(jj-1);
kj_aj = kj*a(jj);

j_kj_aj = sbesselj(n,kj_aj,'norm');
j_kj_aj__ = sbesselj(n,kj_aj__,'norm');
h_kj_aj = sbesselh(n,1,kj_aj,'norm');
h_kj_aj__ = sbesselh(n,1,kj_aj__,'norm');

[R_out_jj T_out_jj] = multilayer_out(mode,u_out_jj,e_out_jj,a_out_jj,k0,n);
[R_std_jj T_std_jj] = multilayer_std(mode,u_std_jj,e_std_jj,a_std_jj,k0,n);

alpha = j_kj_aj__*h_kj_aj/( h_kj_aj__*j_kj_aj ) * exp( 1i*(kj_aj-kj_aj__) + abs(imag(kj_aj__))-abs(imag(kj_aj)) ) ;

%if isnan(alpha) || isinf(alpha)
%	alpha = ( a(jj-1)/a(ii) )^(2*n+1) ;
%end

M_jj = 1/( 1 - R_out_jj(1)*R_std_jj(end) * alpha ); 

%%%%%%%% r and rp are in the same layer %%%%%%%%%%%%%%%%%%%%%
if ii==jj
  %%%%%%%% calculate parameters relevant to r and rp  %%%%%%%%%
  if r>rp
  	rg = r;
  	rl = rp;
  else
  	rg = rp;
  	rl = r;
  end

kj_rl = kj*rl;
kj_rg = kj*rg;

j_kj_rg = sbesselj(n,kj_rg,'norm');
j_kj_rl = sbesselj(n,kj_rl,'norm');
der_j_kj_rg = der_sbesselj(n,kj_rg,'norm');
der_j_kj_rl = der_sbesselj(n,kj_rl,'norm');

h_kj_rg = sbesselh(n,1,kj_rg,'norm');
h_kj_rl = sbesselh(n,1,kj_rl,'norm');
der_h_kj_rg = der_sbesselh(n,1,kj_rg,'norm');
der_h_kj_rl = der_sbesselh(n,1,kj_rl,'norm');


	F = h_kj_rg * j_kj_rl * exp( 1i*kj_rg+abs(imag(kj_rl)) ) ... 
		 + h_kj_rg * R_std_jj(end) * j_kj_aj__ * h_kj_rl / h_kj_aj__ * exp( 1i*kj_rg+1i*kj_rl-1i*kj_aj__+abs(imag(kj_aj__)) ) ...
		 + R_out_jj(1) * h_kj_aj * j_kj_rg / j_kj_aj * j_kj_rl * exp( 1i*kj_aj+abs(imag(kj_rg-kj_aj+kj_rl)) ) ...
		 + R_out_jj(1) * h_kj_aj * j_kj_rg / j_kj_aj * R_std_jj(end) * j_kj_aj__ * h_kj_rl / h_kj_aj__ * exp( 1i*kj_aj+1i*kj_rl-1i*kj_aj__+abs(imag(kj_rg-kj_aj+kj_aj__)) );

	dF_rl = kj * ( h_kj_rg * der_j_kj_rl * exp( 1i*kj_rg+abs(imag(kj_rl)) ) ...
		       + h_kj_rg * R_std_jj(end) * j_kj_aj__ * der_h_kj_rl / h_kj_aj__ * exp( 1i*kj_rg+1i*kj_rl-1i*kj_aj__+abs(imag(kj_aj__)) ) ...
		       + R_out_jj(1) * h_kj_aj * j_kj_rg / j_kj_aj * der_j_kj_rl * exp( 1i*kj_aj+abs(imag(kj_rg-kj_aj+kj_rl)) ) ...
		       + R_out_jj(1) * h_kj_aj * j_kj_rg / j_kj_aj * R_std_jj(end) * j_kj_aj__ * der_h_kj_rl / h_kj_aj__ * exp( 1i*kj_aj+1i*kj_rl-1i*kj_aj__+abs(imag(kj_rg-kj_aj+kj_aj__)) ) );

	dF_rg = kj * ( der_h_kj_rg * j_kj_rl * exp( 1i*kj_rg+abs(imag(kj_rl)) ) ...
		       + der_h_kj_rg * R_std_jj(end) * j_kj_aj__ * h_kj_rl / h_kj_aj__ * exp( 1i*kj_rg+1i*kj_rl-1i*kj_aj__+abs(imag(kj_aj__)) ) ...
		       + R_out_jj(1) * h_kj_aj * der_j_kj_rg / j_kj_aj * j_kj_rl * exp( 1i*kj_aj+abs(imag(kj_rg-kj_aj+kj_rl)) ) ...
		       + R_out_jj(1) * h_kj_aj * der_j_kj_rg / j_kj_aj * R_std_jj(end) * j_kj_aj__ * h_kj_rl / h_kj_aj__ * exp( 1i*kj_aj+1i*kj_rl-1i*kj_aj__+abs(imag(kj_rg-kj_aj+kj_aj__)) ) );

	d2F = kj^2 * ( der_h_kj_rg * der_j_kj_rl * exp( 1i*kj_rg+abs(imag(kj_rl)) ) ...
		       + der_h_kj_rg * R_std_jj(end) * j_kj_aj__ / h_kj_aj__ * der_h_kj_rl * exp( 1i*kj_rg+1i*kj_rl-1i*kj_aj__+abs(imag(kj_aj__)) ) ...
		       + R_out_jj(1) * h_kj_aj / j_kj_aj * der_j_kj_rg * der_j_kj_rl * exp( 1i*kj_aj+abs(imag(kj_rg-kj_aj+kj_rl)) ) ...
		       + R_out_jj(1) * h_kj_aj / j_kj_aj * der_j_kj_rg * R_std_jj(end) * j_kj_aj__ / h_kj_aj__ * der_h_kj_rl * exp( 1i*kj_aj+1i*kj_rl-1i*kj_aj__+abs(imag(kj_rg-kj_aj+kj_aj__)) ) );
  F = F*M_jj;
  dF_rl = dF_rl*M_jj;
  dF_rg = dF_rg*M_jj;
  d2F = d2F*M_jj;
  
  if r>rp
  	dF_rp = dF_rl;
  	dF_r = dF_rg;
  else
  	dF_rp = dF_rg;
  	dF_r = dF_rl;
  end

elseif ii>jj

  [R_out_ii T_out_ii] = multilayer_out(mode,u_out_ii,e_out_ii,a_out_ii,k0,n);

  %%%%%%%%%% parameters relative to ai and ai__ %%%%%%%%%%%%%%
  ki = k0*sqrt( u(ii)*e(ii) );
  ki_ai__ = ki*a(ii-1);
  ki_ai = ki*a(ii);
  j_ki_ai = sbesselj(n,ki_ai,'norm');
  h_ki_ai = sbesselh(n,1,ki_ai,'norm');
  h_ki_ai__ = sbesselh(n,1,ki_ai__,'norm');
  
  %%%%%%%%%% parameters relative to r and rp %%%%%%%%%%%%%%
  kj_rp = kj*rp;
  ki_r = ki*r;

  j_kj_rp = sbesselj(n,kj_rp,'norm');
  der_j_kj_rp = der_sbesselj(n,kj_rp,'norm');
  h_kj_rp = sbesselh(n,1,kj_rp,'norm');
  der_h_kj_rp = der_sbesselh(n,1,kj_rp,'norm');
  
  j_ki_r = sbesselj(n,ki_r,'norm');
  der_j_ki_r = der_sbesselj(n,ki_r,'norm');
  h_ki_r = sbesselh(n,1,ki_r,'norm');
  der_h_ki_r = der_sbesselh(n,1,ki_r,'norm');
  
  %%%%%%%%% calculate Fi  %%%%%%%%%%%%%%%%%
  Fi = h_ki_r/h_ki_ai__ * exp( 1i*ki_r - 1i*ki_ai__ ) + R_out_ii(1) * h_ki_ai / j_ki_ai * j_ki_r / h_ki_ai__ * exp( 1i*ki_ai-abs(imag(ki_ai))+abs(imag(ki_r))-1i*ki*ai__ )  ; 
     
  %%%%%%%%% calculate Fj  %%%%%%%%%%%%%%%%%
  Fj = j_kj_rp*h_kj_aj*exp(abs(imag(kj_rp))+1i*kj_aj) + R_std_jj(end)* j_kj_aj__/h_kj_aj__*h_kj_rp*h_kj_aj*exp(abs(imag(kj_aj__))-1i*kj*aj__+1i*kj_rp+1i*kj_aj) ;
  
  F = Fi*Fj*M_jj*T_out_jj(ii-jj);
  dF_r = NaN;
  dF_rp = NaN;
  d2F = NaN;

else

  [R_std_ii T_std_ii] = multilayer_std(mode,u_std_ii,e_std_ii,a_std_ii,k0,n);

  %%%%%%%%%% parameters relative to ai and ai__ %%%%%%%%%%%%%%
  ki = k0*sqrt( u(ii)*e(ii) );
  ki_ai__ = ki*a(ii-1);
  ki_ai = ki*a(ii);
  j_ki_ai = sbesselj(n,ki_ai,'norm');
  j_ki_ai__ = sbesselj(n,ki_ai__,'norm');
  h_ki_ai__ = sbesselh(n,1,ki_ai__,'norm');

  %%%%%%%%%% parameters relative to r and rp %%%%%%%%%%%%%%
  kj_rp = kj*rp;
  ki_r = ki*r;
  
  j_kj_rp = sbesselj(n,kj_rp,'norm');
  der_j_kj_rp = der_sbesselj(n,kj_rp,'norm');
  h_kj_rp = sbesselh(n,1,kj_rp,'norm');
  der_h_kj_rp = der_sbesselh(n,1,kj_rp,'norm');
  
  j_ki_r = sbesselj(n,ki_r,'norm');
  der_j_ki_r = der_sbesselj(n,ki_r,'norm');
  h_ki_r = sbesselh(n,1,ki_r,'norm');
  der_h_ki_r = der_sbesselh(n,1,ki_r,'norm');

  %%%%%%%%% calculate Fi  %%%%%%%%%%%%%%%%%
  Fi = j_ki_r/j_ki_ai*exp( abs(imag(ki_r))-abs(imag(ki_ai)) ) + R_std_ii(end)* j_ki_ai__/h_ki_ai__*h_ki_r/j_ki_ai*exp( abs(imag(ki_ai__))-1i*ki*ai__+1i*ki_r-abs(imag(ki_ai)) ) ;

  %%%%%%%%% calculate Fj  %%%%%%%%%%%%%%%%%
  Fj = h_kj_rp*j_kj_aj__ * exp( 1i*kj_rp + abs(imag(kj_aj__)) ) + R_out_jj(1) * h_kj_aj / j_kj_aj * j_kj_rp * j_kj_aj__ * exp( 1i*kj_aj-abs(imag(kj_aj))+abs(imag(kj_rp))+abs(imag(kj*aj__)) )  ; 

  F = Fi*Fj*M_jj*T_std_jj(1+ii-jj+end);
  dF_r = NaN;
  dF_rp = NaN;
  d2F = NaN;
     
end

endfunction

