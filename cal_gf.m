% call the green function of spherically layered media ( gf.m/gf.SPrecursion.m )
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

clear all
addpath('./subroutines');
%%%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%
cal = 1; % (1) H_r; (2) E_r.

e = [2 2 1 8 8]; % relative permittivity.
u = [8 8 1 2 2]; % relative permeability.
%e = [1 1 1 1 1]; % relative permittivity.
%u = [1 1 1 1 1]; % relative permeability.
a = [0.1 1 11 20]; % the radius of each interface of neighboring shells.

% cal Hr
r = 0.95;
theta = pi/2-pi/20;
phi = pi/2;
r_p = 1.09:0.01:1.29;
theta_p = pi/2; 
phi_p = pi/2;

%% cal Er
%r = 11.1;
%theta = pi/2;
%phi = pi/2-pi/10;
%r_p = 7.3:0.1:9.3;
%theta_p = pi/2; 
%phi_p = pi/2;

k0 = 1;
precision = 1e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_rp = length(r_p);
field_value = zeros(1,num_rp);

tic
for ii=1:num_rp
 % field_value(ii) = gf(cal,u,e,a,r,theta,phi,r_p(ii),theta_p,phi_p,k0,precision);
  field_value(ii) = gf_SPrecursion(cal,u,e,a,r,theta,phi,r_p(ii),theta_p,phi_p,k0,precision);
end
toc


switch cal
case {1}
	field = 'Hr';
case {2}
	field = 'Er';
endswitch

filename = ['data/' field];
fid = fopen([filename '.txt'],'w');
fprintf(fid,'%f %f %f \n', [r_p; real(field_value); imag(field_value)]);
fclose(fid)

