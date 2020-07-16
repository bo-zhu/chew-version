% calculate the monolayer reflection and transmission of an standing wave
% refer to Prof. Chew's book
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [R_21 T_21] = monolayer_std(mode,u1,u2,e1,e2,a,k0,ALPHA)
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

H_k2a = H(ALPHA,k2a,'none'); 
der_H_k2a = der_H(ALPHA,k2a,'none');

J_k1a = J(ALPHA,k1a,'none');
J_k2a = J(ALPHA,k2a,'none');
der_J_k1a = der_J(ALPHA,k1a,'none');
der_J_k2a = der_J(ALPHA,k2a,'none');

h_k2a = sbesselh(ALPHA,1,k2a,'none');
j_k1a = sbesselj(ALPHA,k1a,'none');
j_k2a = sbesselj(ALPHA,k2a,'none');


	D = sqrt(e1*u2)*J_k1a*der_H_k2a - sqrt(e2*u1)*H_k2a*der_J_k1a;
	R_21 = ( sqrt(e2*u1)*J_k2a*der_J_k1a - sqrt(e1*u2)*J_k1a*der_J_k2a )/D;
	T_21 = 1i*e1*sqrt(u1/e2)/D ;

%if isinf(R_21) || isnan(R_21) || isinf(T_21) || isnan(T_21)
%	R_21 = (e1-e2)/(e1*ALPHA/(ALPHA+1)+e2);
%	T_21 = e1/(e1+e2+e2/ALPHA)*2*(ALPHA+0.5)/ALPHA;
%end

endfunction	
