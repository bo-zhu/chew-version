%% Calculate the static field of a unit charge distribution 
%% at the circle (r_prime,pi/2) in two concentric PEC spherical shells.
%% The outer shell is grounded, while the inner shell is ungrounded.
%% This is the DC field component of a unit point charge moving circularly
%% between the concentric shells.
%% Written by ZHU Bo at Nanjing University. ( email: bzhu@nju.edu.cn )

function [E_r E_theta U] = circle_charge(r1,r2,r,theta,r_prime,precision)
% r1 : radius of the inner shell.
% r2 : radius of the outer shell.
% r : observation distance.
% theta, phi : observation angle.
% r_prime : source distance.
% precision : relative error of two iterations.

addpath('./subroutines');

if r<r1 || r>r2
	printf("the value of 'r' should be r1<r<r2 ! \n");
	return
end
if r_prime<r1 || r_prime>r2
	printf("the value of 'r_prime' should be r1<r_prime<r2 ! \n");
	return
end

%%%%%%%%%%%%%%%%%%%% initialize %%%%%%%%%%%%%%%%
e0 = 8.854e-12;
U = 0;
E_r = 0;
E_theta = 0;
l = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

old_U = 0;
x = cos(theta);
do 
	old_old_U = old_U;
	old_U = U;
	
	switch l
	case {0}
		p0_0 = 1;
		p0_x = 1;
		p0_0_older = 1;
		p0_x_older = 1;
	case {1}
		p0_0 = 0;
		p0_x = x;
		p1_x = -sqrt(1-x^2);
		p0_0_old = 0;
		p0_x_old = x;
		p1_x_older = p1_x;
	case {2}
		p0_0 = legendre_recursion('legendre','degree',0,p0_0_older,p0_0_old,l-1,0);
		p0_x = legendre_recursion('legendre','degree',x,p0_x_older,p0_x_old,l-1,0);
		p1_x = -3*x*sqrt(1-x^2);
		p0_0_older = p0_0_old;
		p0_0_old = p0_0;
		p0_x_older = p0_x_old;
		p0_x_old = p0_x;
		p1_x_old = p1_x;
	otherwise
		p0_0 = legendre_recursion('legendre','degree',0,p0_0_older,p0_0_old,l-1,0);
		p0_x = legendre_recursion('legendre','degree',x,p0_x_older,p0_x_old,l-1,0);
		p1_x = legendre_recursion('legendre','degree',x,p1_x_older,p1_x_old,l-1,1);
		p0_0_older = p0_0_old;
		p0_0_old = p0_0;
		p0_x_older = p0_x_old;
		p0_x_old = p0_x;
		p1_x_older = p1_x_old;
		p1_x_old = p1_x;
	endswitch
		
	f1 = (r_prime/r1)^(2*l+1);
	f2 = (r_prime/r2)^(2*l+1);
	b2 = -l/(l+1)*f1-1 + (l/(l+1)*f2+1)*(1-f1)/(1-f2);
	b1 = (l+0.5)/(l+1)*p0_0/(2*pi*e0*b2);

	% calculate the potential at the observation point (r,theta)
	if r<r_prime
		SUM = b1*(-(r/r1)^l/r1*(r_prime/r1)^l + (r_prime/r)^l/r); 
		U = U + SUM * p0_x;		
		E_r = E_r - ( l*-b1/(r1*r)*(r/r1)^l*(r_prime/r1)^l - (l+1)*b1/r^2*(r_prime/r)^l )*p0_x;
		if l>0
			E_theta = E_theta + SUM*p1_x;
		endif
	else
		SUM = (f1-1)/(1-f2)*b1/r2*(r_prime/r2)^l*(r/r2)^l + (1-f1)/(1-f2)*b1/r*(r_prime/r)^l;
		U = U + SUM * p0_x;		
		E_r = E_r - ( l*(f1-1)/(1-f2)*b1/(r2*r)*(r_prime/r2)^l*(r/r2)^l - (l+1)*(1-f1)/(1-f2)*b1/r^2*(r_prime/r)^l )*p0_x;
		if l>0
			E_theta = E_theta + SUM*p1_x;
		endif
	endif
	
	l = l+1;

until (l>3 && abs( (old_U-U)/U )<precision && abs( (old_old_U-old_U)/old_U )<precision )
printf('Converged after %u iterations!\n', l-1);

% correct U as the inner shell is ungrounded
B0 = 1/(4*pi*e0)*r1/r_prime*(r2-r_prime)/(r1-r2);
U = U + B0*(1/r2-1/r);
E_r = E_r - B0/r^2;
E_theta = E_theta*x/r/sqrt(1-x^2);


endfunction
