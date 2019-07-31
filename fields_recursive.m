% calculate the E and H fields produced by a unit charge moving circularly around the earth in the theta=pi/2 plane.
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

clear all

addpath('./subroutines');

%%%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%
r1 = 0;
r2 = 60;

a = ( 6000+[r1 r2 90 150 450] )*1e3; % the radius of each interface of neighboring shells.
r_prime = 6010e3;
r = 6035e3;
%a = [1 2 3 4 100e3] ; % the radius of each interface of neighboring shells.
%r_prime = 1.5; % the radius of the charge's orbit.
%r = 30e3 ; % r, theta and phi : the observation point. 

theta_prime = pi/2;
theta = pi/2 - 0.01; 
phi = 0;
v = 3; % linear velocity of the charge.
M_truc = 1+3e3; % the truncation frequency = M_truc * Omeg.
cal = 6; % (1) H_r; (2) E_r; (3) H_theta; (4) E_theta; (5) H_phi; (6) E_phi.
precision = 1e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% physical constants %%%%%%%%%%
u0 = 4*pi*1e-7;
e0 = 8.854e-12;
% electron densities of ionosphere layers
ne_F = 1e12;
ne_E = 0.5e11;
ne_D = 1e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Omeg = v/r_prime;
T = 2*pi/Omeg;
m = 0;
cof = 1;
x = cos(theta);
x_prime = cos(theta_prime);
field_value = zeros(M_truc, 1);
field_value_1 = zeros(M_truc, 1);
now = zeros(1,4);
now_prime = zeros(1,4);

idx = find(a>r);
ii = idx(1);
idx = find(a>r_prime);
jj = idx(1);

MM = 1:50:M_truc;
%MM = 51;
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
	
	seq1 = 2*m+1 : 2 : 2*M-1;
	m = M;
	seq2 = seq1 + 1;
	cof = cof * sqrt( prod(seq1./seq2) );

	now(1) = sqrt(2) * cof * (1-x^2)^(M/2);
	now(2) = x * sqrt(2*M+1) * now(1);
%	now(3) = legendre_recursion('sch',x,now(1),now(2),M+1,M);	
%	now(4) = legendre_recursion('sch',x,now(2),now(3),M+2,M);

	now_prime(1) = sqrt(2) * cof * (1-x_prime^2)^(M/2);
	now_prime(2) = x_prime * sqrt(2*M+1) * now_prime(1);
%	now_prime(3) = legendre_recursion('sch',x_prime,now_prime(1),now_prime(2),M+1,M);	
%	now_prime(4) = legendre_recursion('sch',x_prime,now_prime(2),now_prime(3),M+2,M);

	switch cal
	case {1} % calculate H_r

		delta_old = 1e10;
		n = M;
		sp_old_old = now(1);
		sp_old = now(2);
		sp_prime_old_old = now_prime(1);
		sp_prime_old = now_prime(2);
		do	
			[F dF_rp dF_r d2F] = F_tetm('te',u,e,a,r,r_prime,k0,n);

			if n==M	
				Cp_nM = n*x_prime*now_prime(1) / sin(theta_prime);
				delta = (n+0.5) * F * Cp_nM * now(1);	
			elseif n==M+1	
				Cp_nM = ( n*x_prime*now_prime(2) - sqrt(n^2-M^2)*now_prime(1) ) / sin(theta_prime);
				delta = (n+0.5) * F * Cp_nM * now(2);	
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				Cp_nM = ( n*x_prime*sp_prime - sqrt(n^2-M^2)*sp_prime_old ) / sin(theta_prime);
				delta = (n+0.5) * F * Cp_nM * sp;

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value(M) = field_value(M) + delta;
			if ( abs( delta/field_value(M) )<precision && abs( delta_old/field_value(M) )<precision )
				printf('H_r: converged after %d iterations. \n', n-M);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value(M) = field_value(M) * v*kj/(4i*pi*r) * u(jj)/u(ii) * exp(1i*M*phi); 

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
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				delta = (n+0.5) * (F+r_prime*dF_rp) * sp_prime * sp;

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value(M) = field_value(M) + delta;
			if ( abs( delta/field_value(M) )<precision && abs( delta_old/field_value(M) )<precision )
				printf('E_r: converged after %d iterations. \n', n-M+1);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value(M) = field_value(M) * ( 1i*kj ) / ( 4*pi*e(ii)*e0*r ) * exp(1i*M*phi); 

	case {3} % calculate H_theta

		delta_old = 0;
		n = M;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do
			[F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				delta = (n+0.5)/n/(n+1) * (F+r_prime*dF_rp) * now_prime(n-M+1) * now(n-M+1);	
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				delta = (n+0.5)/n/(n+1) * (F+r_prime*dF_rp) * sp_prime * sp;

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value(M) = field_value(M) + delta;
			if ( abs( delta/field_value(M) )<precision && abs( delta_old/field_value(M) )<precision )
				printf('H_theta_term1: converged after %d iterations. \n', n-M+1);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value(M) = field_value(M) * 1i*M/sin(theta) * w/(4*pi)*kj * exp(1i*M*phi); 

		delta_old = 0;
		n = M+1;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do
			[F dF_rp dF_r d2F] = F_tetm('te',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				C_nM = ( n*x*now(n-M+1) - sqrt(n^2-M^2)*now(n-M) ) / sin(theta);
				delta = (n+0.5)*sqrt(n^2-M^2)/n/(n+1) * (F+r*dF_r) * now_prime(n-M) * C_nM;	
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				C_nM = ( n*x*sp - sqrt(n^2-M^2)*sp_old ) / sin(theta);
				delta = (n+0.5)*sqrt(n^2-M^2)/n/(n+1) * (F+r*dF_r) * sp_prime_old * C_nM;	

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value_1(M) = field_value_1(M) + delta;
			if ( abs( delta/field_value_1(M) )<precision && abs( delta_old/field_value_1(M) )<precision )
				printf('H_theta_term2: converged after %d iterations. \n', n-M+2);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		
		field_value_1(M) = field_value_1(M) * ( 1i*v*kj ) / ( 4*pi*r ) * u(jj)/u(ii) * exp(1i*M*phi); 
		field_value(M) = field_value(M) + field_value_1(M);

	case {4} % calculate E_theta

		delta_old =0;
		n = M+1;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do	
			[F dF_rp dF_r d2F] = F_tetm('te',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				delta = -(n+0.5)*sqrt(n^2-M^2)/n/(n+1) * F * now_prime(n-M) * now(n-M+1);	
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				delta = -(n+0.5)*sqrt(n^2-M^2)/n/(n+1) * F * sp_prime_old * sp;

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value(M) = field_value(M) + delta;
			if ( abs( delta/field_value(M) )<precision && abs( delta_old/field_value(M) )<precision )
				printf('E_theta_term1: converged after %d iterations. \n', n-M+2);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value(M) = field_value(M) * 1i*M / sin(theta)  * w*v*u(jj)*u0*kj/(4*pi) * exp(1i*M*phi) ;

		delta_old = 0;
		n = M;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do
			[F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				if n == M
					C_nM = ( n*x*now(n-M+1) ) / sin(theta);	
				else
					C_nM = ( n*x*now(n-M+1) - sqrt(n^2-M^2)*now(n-M) ) / sin(theta);
				end
				delta = (n+0.5)/n/(n+1) * (F+r*dF_r+r_prime*dF_rp+r*r_prime*d2F) *...
					now_prime(n-M+1) * C_nM	;
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				C_nM = ( n*x*sp - sqrt(n^2-M^2)*sp_old ) / sin(theta);
				delta = (n+0.5)/n/(n+1) * (F+r*dF_r+r_prime*dF_rp+r*r_prime*d2F) *...
					sp_prime * C_nM	;

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value_1(M) = field_value_1(M) + delta;
			if ( abs( delta/field_value_1(M) )<precision && abs( delta_old/field_value_1(M) )<precision )
				printf('E_theta_term2: converged after %d iterations. \n', n-M+1);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value_1(M) = field_value_1(M) * 1i*kj / ( 4*pi*r*e(ii)*e0 ) * exp(1i*M*phi); 
		field_value(M) = field_value(M) + field_value_1(M); 

	case {5} % calculate H_phi

		delta_old = 0;
		n = M;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do
			[F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				if n == M
					C_nM = n*x*now(n-M+1) / sin(theta);	
				else
					C_nM = ( n*x*now(n-M+1) - sqrt(n^2-M^2)*now(n-M) ) / sin(theta);
				end
				delta = (n+0.5)/n/(n+1) * (F+r_prime*dF_rp) * now_prime(n-M+1) * C_nM;	
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				C_nM = ( n*x*sp - sqrt(n^2-M^2)*sp_old ) / sin(theta);
				delta = (n+0.5)/n/(n+1) * (F+r_prime*dF_rp)* sp_prime * C_nM;	

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value(M) = field_value(M) + delta;
			if ( abs( delta/field_value(M) )<precision && abs( delta_old/field_value )<precision )
				printf('H_phi_term1: converged after %d iterations. \n', n-M+1);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value(M) = -field_value(M) * w*kj / ( 4*pi ) * exp(1i*M*phi); 

		delta_old = 0;
		n = M+1;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do
			[F dF_rp dF_r d2F] = F_tetm('te',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				delta = (n+0.5)*sqrt(n^2-M^2)/n/(n+1) * (F+r*dF_r) * now_prime(n-M) * now(n-M+1);	
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				delta = (n+0.5)*sqrt(n^2-M^2)/n/(n+1) * (F+r*dF_r) * sp_prime_old * sp;	

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value_1(M) = field_value_1(M) + delta;
			if ( abs( delta/field_value_1(M) )<precision && abs( delta_old/field_value_1(M) )<precision )
				printf('H_phi_term2: converged after %d iterations. \n', n-M+2);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value_1(M) = field_value_1(M) * ( M*v*kj ) / ( 4*pi*r*sin(theta) ) * u(jj)/u(ii) * exp(1i*M*phi); 
		field_value(M) = field_value(M) - field_value_1(M);

	case {6} % calculate E_phi

		delta_old = 0;
		n = M+1;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do
			[F dF_rp dF_r d2F] = F_tetm('te',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				C_nM = ( n*x*now(n-M+1) - sqrt(n^2-M^2)*now(n-M) ) / sin(theta);
				delta = (n+0.5)*sqrt(n^2-M^2)/n/(n+1) * F * now_prime(n-M) * C_nM;	
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				C_nM = ( n*x*sp - sqrt(n^2-M^2)*sp_old ) / sin(theta) ;
				delta = (n+0.5)*sqrt(n^2-M^2)/n/(n+1) * F * sp_prime_old * C_nM;	

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value(M) = field_value(M) + delta;
			if ( abs( delta/field_value(M) )<precision && abs( delta_old/field_value(M) )<precision )
				printf('E_phi_term1: converged after %d iterations. \n', n-M+2);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value(M) = field_value(M) * ( w*v*u(jj)*u0*kj ) / ( 4*pi ) * exp(1i*M*phi); 

		delta_old = 0;
		n = M;
		sp_old_old = now(3);
		sp_old = now(4);
		sp_prime_old_old = now_prime(3);
		sp_prime_old = now_prime(4);
		do
			[F dF_rp dF_r d2F] = F_tetm('tm',u,e,a,r,r_prime,k0,n);

			if n<=M+3	
				delta = (n+0.5)/n/(n+1) * (F + r*dF_r + r_prime*dF_rp + r*r_prime*d2F) *...
					 now_prime(n-M+1) * now(n-M+1);	
			else
				sp = legendre_recursion('sch',x,sp_old_old,sp_old,n-1,M); 
				sp_prime = legendre_recursion('sch',x_prime,sp_prime_old_old,sp_prime_old,n-1,M); 
				delta = (n+0.5)/n/(n+1) * (F + r*dF_r + r_prime*dF_rp + r*r_prime*d2F) *...
					 sp_prime * sp;	

				sp_old_old = sp_old;
				sp_old = sp;
				sp_prime_old_old = sp_prime_old;
				sp_prime_old = sp_prime;
			end

			field_value_1(M) = field_value_1(M) + delta;
			if ( abs( delta/field_value_1(M) )<precision && abs( delta_old/field_value_1(M) )<precision )
				printf('E_phi_term2: converged after %d iterations. \n', n-M+1);
				break
			else
				n = n+1; 
				delta_old = delta ;
			endif
		until (0)
		field_value_1(M) = field_value_1(M) * M*kj*exp(1i*M*phi) / ( 4*pi*e(ii)*e0*r*sin(theta) ) ; 
		field_value(M) = field_value(M) - field_value_1(M); 

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


