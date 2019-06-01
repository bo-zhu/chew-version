% calculate the E and H fields produced by a unit charge moving circularly around the earth in the theta=pi/2 plane.
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

clear all

addpath('./subroutines');

%%%%%%%% physical constants %%%%%%%%%%
u0 = 4*pi*1e-7;
e0 = 8.854e-12;
% electron densities of ionosphere layers
ne_F = 1e12;
ne_E = 0.5e11;
ne_D = 1e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%
r1 = 0;
r2 = 60;
a = ( 6000+[r1 r2 90 150 450] )*1e3; % the radius of each interface of neighboring shells.
r_prime = 6000e3 + 10e3; % the radius of the charge's orbit.
r = 6000e3 + 35e3; % r, theta and phi : the observation point. 
theta_prime = pi/2;
theta = pi/2;
phi = 0;
v = 300; % linear velocity of the charge.
M_truc = 1+3e3; % the truncation frequency = M_truc * Omeg.
cal = 2; % (1) H_r; (2) E_r; (3) H_theta; (4) E_theta; (5) H_phi; (6) E_phi.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = find(a>r);
ii = idx(1);
idx = find(a>r_prime);
jj = idx(1);

%%%%%%%%%% initial sp(x) %%%%%%%%%%%
x = cos(theta);
sp0_2 =(3*x^2-1)/sqrt(2);
sp0_3 = x*(5*x^2-3)/sqrt(2);
sp1_1 = sqrt(1-x^2);
sp1_2 = sqrt(3)*x*sqrt(1-x^2); 
old = [0 0 sp0_2 sp0_3];
now = [sp1_1 sp1_2 0 0];
now(3) = legendre_recursion('sch','degree',x,now(1),now(2),2,1);	
now(4) = legendre_recursion('sch','degree',x,now(2),now(3),3,1);
new = zeros(1,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% initial sp(x_prime) %%%%%%%%
x_prime = cos(theta_prime);
sp0_2_prime =(3*x_prime^2-1)/sqrt(2);
sp0_3_prime = x_prime*(5*x_prime^2-3)/sqrt(2);
sp1_1_prime = sqrt(1-x_prime^2);
sp1_2_prime = sqrt(3)*x_prime*sqrt(1-x_prime^2); 
old_prime = [0 0 sp0_2_prime sp0_3_prime];
now_prime = [sp1_1_prime sp1_2_prime 0 0];
now_prime(3) = legendre_recursion('sch','degree',x,now_prime(1),now_prime(2),2,1);	
now_prime(4) = legendre_recursion('sch','degree',x,now_prime(2),now_prime(3),3,1);
new_prime = zeros(1,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Omeg = v/r_prime;
T = 2*pi/Omeg;
m = 1;
MM = 1:50:M_truc;
field_value = zeros(M_truc, 1);

for M = MM 
	w = M*Omeg;

	%%%%%%%% media distribution 1 %%%%%%%%%%
	% the order of the elements in u and e shall be from +z to -z direction 
	u = [1 1 1 1 1 1]; 
	e_F = plasma_para(ne_F,w);
	e_E = plasma_para(ne_E,w);
	e_D = plasma_para(ne_D,w);
	e_earth = 1 + 1i*1e9/e0/w;
	e = [e_earth 1 e_D e_E e_F 1];
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	%%%%%%%% media distribution 2 %%%%%%%%%%
%	% the order of the elements in u and e shall be from +z to -z direction 
%	u = [1 1 1 1 1 1]; 
%	e_earth = 1 + 1i*1e9/e0/w;
%	e = [e_earth 1 1 1 1 1];
%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	%%%%%%%% media distribution 3 %%%%%%%%%%
%	% the order of the elements in u and e shall be from +z to -z direction 
%	u = [1 1 1 1 1 1]; 
%	e = [1 1 1 1 1 1];
%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	k0 = w*sqrt(e0*u0);
	kj = k0*sqrt(e(jj)*u(jj));

	while (m<M)
		new(1) = legendre_recursion('sch','order',x,old(3),now(2),m+1,m);
		new(2) = legendre_recursion('sch','order',x,old(4),now(3),m+2,m);
		new(3) = legendre_recursion('sch','degree',x,new(1),new(2),m+2,m+1);	
		new(4) = legendre_recursion('sch','degree',x,new(2),new(3),m+3,m+1);

		new_prime(1) = legendre_recursion('sch','order',x,old_prime(3),now_prime(2),m+1,m);
		new_prime(2) = legendre_recursion('sch','order',x,old_prime(4),now_prime(3),m+2,m);
		new_prime(3) = legendre_recursion('sch','degree',x,new_prime(1),new_prime(2),m+2,m+1);	
		new_prime(4) = legendre_recursion('sch','degree',x,new_prime(2),new_prime(3),m+3,m+1);

		m = m+1;
		old = now;
		now = new;
		old_prime = now_prime;
		now_prime = new_prime;
	endwhile

	switch cal
	case {1} % calculate H_r

		delta_old =0;
		n = M+1;
		do	
			[F dF_rp dF_r d2F] = F_tetm('te',u,e,a,r,r_prime,k0,n);

			switch n	
			case {M+1}
				delta = sqrt(n^2-M^2) * F * now_prime(1) * now(2);	
				sp_old_old = now(2);
				sp_prime_old_old = now_prime(1);
			case {M+2}
				delta = sqrt(n^2-M^2) * F * now_prime(2) * now(3);	
				sp_old = now(2);
				sp_prime_old = now_prime(2);
			otherwise
				sp = legendre_recursion('sch','degree',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch','degree',x,sp_prime_old_old,sp_prime_old,n-1,M); 
				delta = (n+0.5) * (F+r_prime*dF_rp) * sp_prime * sp;

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			endswitch

			field_value(M) = field_value(M) + delta;
			if ( abs( delta/field_value(M) )<1e-4 && abs( delta_old/field_value(M) )<1e-4 )
				printf('TE: converged after %d iterations. \n', n-M+1);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value(M) = field_value(M) * ( 1i*v*kj*u(jj) ) / ( 4*pi*r*u(ii) ) * exp(1i*M*phi); 

	case {2} % calculate E_r

		delta_old = 0;
		n = M;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do
			[F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				delta = (n+0.5) * (F+r_prime*dF_rp) * now_prime(n-M+1) * now(n-M+1);	
			else
				sp = legendre_recursion('sch','degree',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch','degree',x,sp_prime_old_old,sp_prime_old,n-1,M); 
				delta = (n+0.5) * (F+r_prime*dF_rp) * sp_prime * sp;

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			%switch n	
			%case {M}
			%	delta = (n+0.5) * (F+r_prime*dF_rp) * now_prime(1) * now(1);	
			%	sp_old_old = now(1);
			%	sp_prime_old_old = now_prime(1);
			%case {M+1}
			%	delta = (n+0.5) * (F+r_prime*dF_rp) * now_prime(2) * now(2);	
			%	sp_old = now(2);
			%	sp_prime_old = now_prime(2);
			%otherwise
			%	sp = legendre_recursion('sch','degree',x,sp_old_old,sp_old,n-1,M); 
			%	sp_prime = legendre_recursion('sch','degree',x,sp_prime_old_old,sp_prime_old,n-1,M); 
			%	delta = (n+0.5) * (F+r_prime*dF_rp) * sp_prime * sp;

			%	sp_old_old = sp_old;
			%	sp_old = sp;
			%	sp_prime_old_old = sp_prime_old;
			%	sp_prime_old = sp_prime;
			%endswitch

			field_value(M) = field_value(M) + delta;
			if ( abs( delta/field_value(M) )<1e-4 && abs( delta_old/field_value(M) )<1e-4 )
				printf('TM: converged after %d iterations. \n', n-M+1);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value(M) = field_value(M) * ( 1i*kj ) / ( 4*pi*e(ii)*e0*r ) * exp(1i*M*phi); 

	case {3} % calculate H_theta
	case {4} % calculate E_theta
	case {5} % calculate H_phi
	case {6} % calculate E_phi
	endswitch
endfor


switch cal
case {1}
	field = 'Hr';
case {2}
	field = 'Er';
case {3}
	field = 'Htheta';
case {4}
	field = 'Etheta';
case {5}
	field = 'Hphi';
case {6}
	field = 'Ephi';
endswitch

filename = ['data/' field '_10km_' num2str((r-6000e3)/1e3) 'km_' num2str(v) 'mps'];
fid = fopen([filename '_w.txt'],'w');
fprintf(fid,'%10.5f %10.5f \n', [real(field_value(MM))'; imag(field_value(MM))']);
fclose(fid)

field_value = interp1(MM,field_value(MM),(1:M_truc)');

dt = 0.1*T/M_truc;
t = 0:dt:100*dt;
len_t = length(t);
M = repmat( (1:M_truc)', 1, len_t);
t = repmat(t, M_truc, 1);
wt = (M*Omeg).*t;

field_value_abs = repmat(abs(field_value), 1, len_t);
field_value_angle = repmat(angle(field_value), 1, len_t);
field_value_t = 2*sum( field_value_abs.*cos(field_value_angle-wt) , 1 );

figure(1)
plot(t(1,:),field_value_t);
xlabel('t (s)');
ylabel(field);
saveas(1,[filename '_t.pdf']);

fid = fopen([filename '_t.txt'],'w');
fprintf(fid,'%10.5f %10.5f\n', [t(1,:) ; field_value_t]);
fclose(fid)


