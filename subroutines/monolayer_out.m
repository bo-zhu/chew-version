% calculate the monolayer reflection and transmission of an outging wave
% refer to Prof. Chew's book
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [R_12 T_12] = monolayer_out(mode,u1,u2,e1,e2,a,k0,ALPHA)
% mode : TE or TM
% u,e : relative parameters. 1(2) denotes the inner(outer) layer.
% a : the radius of the inner sphere
% k0 : free sapce wave number, w*sqrt(e0*u0).
% ALPHA : degree of the Hankel function.
% exp(-iwt) is used.


k1a = k0*sqrt(u1*e1)*a;
k2a = k0*sqrt(u2*e2)*a;

e0 = 8.854e-12;
u0 = 4*pi*1e-7;
e1 = e0*e1;
e2 = e0*e2;
u1 = u0*u1;
u2 = u0*u2;

if mode=='TE' || mode=='te'  % duality principle
	temp = u1;
	u1 = e1;
	e1 = temp;
	temp = u2;
	u2 = e2;
	e2 = temp;
endif

H_k1a = H(ALPHA,k1a,'none'); 
H_k2a = H(ALPHA,k2a,'none'); 
der_H_k1a = der_H(ALPHA,k1a,'none');
der_H_k2a = der_H(ALPHA,k2a,'none');

J_k1a = J(ALPHA,k1a,'none');
der_J_k1a = der_J(ALPHA,k1a,'none');

h_k1a = sbesselh(ALPHA,1,k1a,'none');
h_k2a = sbesselh(ALPHA,1,k2a,'none');
j_k1a = sbesselj(ALPHA,k1a,'none');

	D = sqrt(e1*u2)*J_k1a*der_H_k2a - sqrt(e2*u1)*H_k2a*der_J_k1a;
	R_12 = ( sqrt(e2*u1)*H_k2a*der_H_k1a - sqrt(e1*u2)*H_k1a*der_H_k2a )/D;
	T_12 = 1i*e2*sqrt(u2/e1)/D ;

%if isinf(R_12) || isnan(R_12) || isinf(T_12) || isnan(T_12)
%	R_12 = (e2-e1)/(e1+e2+e2/ALPHA);
%	T_12 = e2/(e1+e2+e2/ALPHA)*2*(ALPHA+0.5)/ALPHA;
%end

endfunction	
