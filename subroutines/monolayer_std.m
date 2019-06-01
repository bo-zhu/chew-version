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

H_k2a = H(ALPHA,k2a,'norm'); 
der_H_k2a = der_H(ALPHA,k2a,'norm');

J_k1a = J(ALPHA,k1a,'norm');
J_k2a = J(ALPHA,k2a,'norm');
der_J_k1a = der_J(ALPHA,k1a,'norm');
der_J_k2a = der_J(ALPHA,k2a,'norm');

h_k2a = sbesselh(ALPHA,1,k2a,'norm');
j_k1a = sbesselj(ALPHA,k1a,'norm');
j_k2a = sbesselj(ALPHA,k2a,'norm');

if abs(k1a)<1 
	if abs(k2a)<1
		R_21 = (e1-e2)/(e1*ALPHA/(ALPHA+1)+e2);
		T_21 = e1/(e1+e2+e2/ALPHA)*2*(ALPHA+0.5)/ALPHA;
	else
		D = sqrt(e1*u2)*k1a*der_H_k2a - sqrt(e2*u1)*(ALPHA+1)*H_k2a;
		R_21 = h_k2a/j_k2a*( sqrt(e2*u1)*(ALPHA+1)*J_k2a - sqrt(e1*u2)*k1a*der_J_k2a )/D;
		T_21 = 1i*e1*sqrt(u1/e2)/j_k2a/D * exp(-1i*k2a-abs(imag(k2a)));
	endif
else
	if abs(k2a)<1
		D = sqrt(e1*u2)*ALPHA*J_k1a + sqrt(e2*u1)*k2a*der_J_k1a;
		R_21 = ( sqrt(e1*u2)*(ALPHA+1)*J_k1a - sqrt(e2*u1)*k2a*der_J_k1a )/D;
		T_21 = e1*sqrt(u1/e2)*(2*ALPHA+1)*k2a*j_k1a/D;
	else 
		D = sqrt(e1*u2)*J_k1a*der_H_k2a - sqrt(e2*u1)*H_k2a*der_J_k1a;
		R_21 = h_k2a/j_k2a*( sqrt(e2*u1)*J_k2a*der_J_k1a - sqrt(e1*u2)*J_k1a*der_J_k2a )/D;
		T_21 = j_k1a/j_k2a*1i*e1*sqrt(u1/e2)/D * exp(-1i*k2a-abs(imag(k2a)));
	endif
endif	

endfunction	
