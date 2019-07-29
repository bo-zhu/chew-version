%% Calculate the static field of a unit charge distribution 
%% at the circle (r_prime,theta_prime) in two concentric PEC spherical shells.
%% The outer shell is grounded, while the inner shell is ungrounded.
%% This is the DC field component of a unit point charge moving circularly
%% between the concentric shells.
%% Written by ZHU Bo at Nanjing University. ( email: bzhu@nju.edu.cn )

function [E_r E_theta U] = circle_charge(r1,r2,r,theta,r_prime,theta_prime,precision)
% r1 : radius of the inner shell.
% r2 : radius of the outer shell.
% r : observation distance.
% theta : observation angle.
% r_prime : source distance.
% theta_prime: source angle.
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
l = 0;

U = 0;
old_U = 0;
U_converge = 0;
U_zero = 0;

E_r = 0;
old_E_r = 0;
E_r_converge = 0;
E_r_zero = 0;

E_theta = 0;
old_E_theta = 0;
E_theta_converge = 0;
E_theta_zero = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = cos(theta);
xp = cos(theta_prime);
do 
	
	switch l
	case {0}
		p0_xp = 1;
		p0_x = 1;
		p0_xp_older = 1;
		p0_x_older = 1;
	case {1}
		p0_xp = xp;
		p0_x = x;
		p1_x = -sqrt(1-x^2);
		p0_xp_old = xp;
		p0_x_old = x;
		p1_x_older = p1_x;
	case {2}
		p0_xp = legendre_recursion('legendre','degree',xp,p0_xp_older,p0_xp_old,l-1,0);
		p0_x = legendre_recursion('legendre','degree',x,p0_x_older,p0_x_old,l-1,0);
		p1_x = -3*x*sqrt(1-x^2);
		p0_xp_older = p0_xp_old;
		p0_xp_old = p0_xp;
		p0_x_older = p0_x_old;
		p0_x_old = p0_x;
		p1_x_old = p1_x;
	otherwise
		p0_xp = legendre_recursion('legendre','degree',0,p0_xp_older,p0_xp_old,l-1,0);
		p0_x = legendre_recursion('legendre','degree',x,p0_x_older,p0_x_old,l-1,0);
		p1_x = legendre_recursion('legendre','degree',x,p1_x_older,p1_x_old,l-1,1);
		p0_xp_older = p0_xp_old;
		p0_xp_old = p0_xp;
		p0_x_older = p0_x_old;
		p0_x_old = p0_x;
		p1_x_older = p1_x_old;
		p1_x_old = p1_x;
	endswitch
		
	f1 = (r_prime/r1)^(2*l+1);
	f2 = (r_prime/r2)^(2*l+1);
	b2 = -l/(l+1)*f1-1 + (l/(l+1)*f2+1)*(1-f1)/(1-f2);
	b1 = (l+0.5)/(l+1)*p0_xp/(2*pi*e0*b2);

	% calculate the potential at the observation point (r,theta)
	switch r<r_prime
        case 1
		SUM = b1*(-(r/r1)^l/r1*(r_prime/r1)^l + (r_prime/r)^l/r); 
		
		if U_converge == 0
			old_old_U = old_U;
			old_U = U;
			U = U + SUM * p0_x;		
			if (l>3 && abs( (old_U-U)/U )<precision && abs( (old_old_U-old_U)/old_U )<precision )
				U_converge = 1;
				printf(' U is converged after %u iterations!\n', l-1);
			end
			if (l>3 && U == 0 & old_U == 0)
				U_zero = U_zero + 1;
				if U_zero > iter_max
					U_converge = 1;
					printf(' U seems to be ZERO.\n ');
				end
			end
		end
		
		if E_r_converge == 0
			old_old_E_r = old_E_r;
			old_E_r = E_r;
			E_r = E_r - ( l*-b1/(r1*r)*(r/r1)^l*(r_prime/r1)^l - (l+1)*b1/r^2*(r_prime/r)^l )*p0_x;
			if (l>3 && abs( (old_E_r-E_r)/E_r )<precision ...
				&& abs( (old_old_E_r-old_E_r)/old_E_r )<precision )
				E_r_converge = 1;
				printf(' E_r is converged after %u iterations!\n', l-1);
			end
			if (l>3 && E_r == 0 & old_E_r == 0)
				E_r_zero = E_r_zero + 1;
				if E_r_zero > iter_max
					E_r_converge = 1;
					printf(' E_r seems to be ZERO. \n');
				end
			end
		end

		if ( l>0 && E_theta_converge == 0 )
			old_old_E_theta = old_E_theta;
			old_E_theta = E_theta;
			E_theta = E_theta + SUM*p1_x;
			if (l>3 && abs( (old_E_theta-E_theta)/E_theta )<precision ...
				&& abs( (old_old_E_theta-old_E_theta)/old_E_theta )<precision )
				E_theta_converge = 1;
				printf(' E_theta is converged after %u iterations!\n', l-1);
			end
			if (l>3 && E_theta == 0 & old_E_theta == 0)
				E_theta_zero = E_theta_zero + 1;
				if E_theta_zero > iter_max
					E_theta_converge = 1;
					printf(' E_theta seems to be ZERO. \n');
				end
			end
		end

	case 0
		SUM = (f1-1)/(1-f2)*b1/r2*(r_prime/r2)^l*(r/r2)^l + (1-f1)/(1-f2)*b1/r*(r_prime/r)^l;
		
		if U_converge == 0
			old_old_U = old_U;
			old_U = U;
keyboard
			U = U + SUM * p0_x;		
			if (l>3 && abs( (old_U-U)/U )<precision && abs( (old_old_U-old_U)/old_U )<precision )
				U_converge = 1;
				printf(' U is converged after %u iterations!\n', l-1);
			end
			if (l>3 && U == 0 & old_U == 0)
				U_zero = U_zero + 1;
				if U_zero > iter_max
					U_converge = 1;
					printf(' U seems to be ZERO.\n ');
				end
			end
		end
		
		if E_r_converge == 0
			old_old_E_r = old_E_r;
			old_E_r = E_r;
			E_r = E_r - ( l*(f1-1)/(1-f2)*b1/(r2*r)*(r_prime/r2)^l*(r/r2)^l - ...
				      (l+1)*(1-f1)/(1-f2)*b1/r^2*(r_prime/r)^l )*p0_x;
			if (l>3 && abs( (old_E_r-E_r)/E_r )<precision ...
				&& abs( (old_old_E_r-old_E_r)/old_E_r )<precision )
				E_r_converge = 1;
				printf(' E_r is converged after %u iterations!\n', l-1);
			end
			if (l>3 && E_r == 0 & old_E_r == 0)
				E_r_zero = E_r_zero + 1;
				if E_r_zero > iter_max
					E_r_converge = 1;
					printf(' E_r seems to be ZERO. \n');
				end
			end
		end
	
		if ( l>0 && E_theta_converge == 0 )
			old_old_E_theta = old_E_theta;
			old_E_theta = E_theta;
			E_theta = E_theta + SUM*p1_x;
			if (l>3 && abs( (old_E_theta-E_theta)/E_theta )<precision ...
				&& abs( (old_old_E_theta-old_E_theta)/old_E_theta )<precision )
				E_theta_converge = 1;
				printf(' E_theta is converged after %u iterations!\n', l-1);
			end
			if (l>3 && E_theta == 0 & old_E_theta == 0)
				E_theta_zero = E_theta_zero + 1;
				if E_theta_zero > iter_max
					E_theta_converge = 1;
					printf(' E_theta seems to be ZERO. \n');
				end
			end
		end
	end
	
	l = l+1;

until (U_converge && E_theta_converge && E_r_converge)

% correct U as the inner shell is ungrounded
B0 = 1/(4*pi*e0)*r1/r_prime*(r2-r_prime)/(r1-r2);
U = U + B0*(1/r2-1/r);
E_r = E_r - B0/r^2;
E_theta = E_theta*x/r/sqrt(1-x^2);


endfunction
