% calculate the multilayer reflection and transmission of an standing wave
% refer to Prof. Chew's book
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [R_std T_std] = multilayer_std(mode,u,e,a,k0,ALPHA)
% mode : TE or TM
% u,e : relative parameters of each layer. From inner to outer.
% a : the radius of each interface. From inner to outer.
% k0 : free sapce wave number, w*sqrt(e0*u0).
% ALPHA : degree of the Hankel function.
% R,T : reflection and transmission of each interface. Stored from inner to outer.
% exp(-iwt) is used.


N = length(a) + 1; % the total number of regions, equal to the lengths of u and e.
R_std = zeros(1,N-1);
T_std = zeros(1,N-1);
S = zeros(1,N-1);
jn = zeros(1,N-2);

[R_std(1) S(1)] = monolayer_std(mode,u(1),u(2),e(1),e(2),a(1),k0,ALPHA);
if N>2
	for ii=2:1:N-1
		[R_12 T_12] = monolayer_out(mode,u(ii),u(ii+1),e(ii),e(ii+1),a(ii),k0,ALPHA);
		[R_21 T_21] = monolayer_std(mode,u(ii),u(ii+1),e(ii),e(ii+1),a(ii),k0,ALPHA);
		k1 = k0*sqrt( u(ii)*e(ii) );
		k1a1 = k1*a(ii);
		k1a0 = k1*a(ii-1);

		h_k1a1 = sbesselh(ALPHA,1,k1a1,'norm');
		h_k1a0 = sbesselh(ALPHA,1,k1a0,'norm');
		j_k1a1 = sbesselj(ALPHA,k1a1,'norm');
		j_k1a0 = sbesselj(ALPHA,k1a0,'norm');

		if abs(k1a1)<1
			alpha = (a(ii-1)/a(ii))^(2*ALPHA+1);
			jn(ii-1) = (k1a0/k1a1)^ALPHA;
		elseif abs(k1a0)<1
			Factor = 2*ALPHA-1:-2:1;
			PROD = prod(k1a0./Factor);
			alpha = 1i*pi*k1a0/(2*ALPHA+1) * ( PROD/gamma(0.5) )^2 * h_k1a1/j_k1a1 * exp(1i*k1a1-abs(imag(k1a1)));
			jn(ii-1) = j_k1a0/j_k1a1 * exp( abs(imag(k1a0))-abs(imag(k1a1)) );
		else
			alpha = ( h_k1a1*j_k1a0 )/( h_k1a0*j_k1a1 ) * exp(1i*k1a1-1i*k1a0+abs(imag(k1a0))-abs(imag(k1a1))); 
			jn(ii-1) = j_k1a0/j_k1a1 * exp( abs(imag(k1a0))-abs(imag(k1a1)) );
		end
		if isnan(alpha) || isinf(alpha) || isnan(jn(ii-1)) || isinf(jn(ii-1))
			alpha = (a(ii)/a(ii+1))^(2*ALPHA+1);
			jn(ii-1) = (k1a0/k1a1)^ALPHA;
		end
		D = 1-R_12*R_std(ii-1)*alpha;
		R_std(ii) = R_21 + T_12*T_21*R_std(ii-1)*alpha/D;
		S(ii) = T_21/D;
	end
	T_std(1:N-2) = fliplr(cumprod( fliplr(jn.*S(2:N-1)) )).*S(1:N-2);
endif

T_std(N-1) = S(N-1);



endfunction	
