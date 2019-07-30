% calculate the multilayer reflection and transmission of an outging wave
% refer to Prof. Chew's book
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [R_out T_out] = multilayer_out(mode,u,e,a,k0,ALPHA)
% mode : TE or TM
% u,e : relative parameters of each layer. From inner to outer.
% a : the radius of each interface. From inner to outer.
% k0 : free sapce wave number, w*sqrt(e0*u0).
% ALPHA : degree of the Hankel function.
% R,T : reflection and transmission of each interface. Stored from inner to outer.
% exp(-iwt) is used.


N = length(a) + 1; % the total number of regions, equal to the lengths of u and e.
R_out = zeros(1,N-1);
T_out = zeros(1,N-1);
S = zeros(1,N-1);
hn = zeros(1,N-2);

[R_out(N-1) S(N-1)] = monolayer_out(mode,u(N-1),u(N),e(N-1),e(N),a(N-1),k0,ALPHA);
if N>2
	for ii=N-2:-1:1
		[R_12 T_12] = monolayer_out(mode,u(ii),u(ii+1),e(ii),e(ii+1),a(ii),k0,ALPHA);
		[R_21 T_21] = monolayer_std(mode,u(ii),u(ii+1),e(ii),e(ii+1),a(ii),k0,ALPHA);
		k2 = k0*sqrt( u(ii+1)*e(ii+1) ); 
		k2a1 = k2*a(ii);
		k2a2 = k2*a(ii+1);

		h_k2a1 = sbesselh(ALPHA,1,k2a1,'norm');
		h_k2a2 = sbesselh(ALPHA,1,k2a2,'norm');
		j_k2a1 = sbesselj(ALPHA,k2a1,'norm');
		j_k2a2 = sbesselj(ALPHA,k2a2,'norm');

		if abs(k2a2)<1
			alpha = (a(ii)/a(ii+1))^(2*ALPHA+1);
			hn(ii) = (k2a2/k2a1)^(-ALPHA-1);
		elseif abs(k2a1)<1
			Factor = 2*ALPHA-1:-2:1;
			SUM = sum(log(k2a1./Factor),'extra');
			alpha = 1i*pi*k2a1/(2*ALPHA+1) * ( exp(SUM)/gamma(0.5) )^2 * h_k2a2/j_k2a2 * exp(1i*k2a2-abs(imag(k2a2)));
			hn(ii) = 1i*sqrt(pi)/gamma(0.5)*h_k2a2*k2a1*exp(1i*k2a2+SUM);
		else
			alpha = ( h_k2a2*j_k2a1 )/( h_k2a1*j_k2a2 ) * exp(1i*k2a2-1i*k2a1+abs(imag(k2a1))-abs(imag(k2a2))); 
			hn(ii) = h_k2a2/h_k2a1 * exp(1i*k2a2-1i*k2a1);
		end
		
		if isnan(alpha) || isinf(alpha) || isnan(hn(ii)) || isinf(hn(ii))
			alpha = (a(ii)/a(ii+1))^(2*ALPHA+1);
			hn(ii) = (k2a2/k2a1)^(-ALPHA-1);
		end

		D = 1-R_21*R_out(ii+1)*alpha;
		R_out(ii) = R_12 + T_21*T_12*R_out(ii+1)*alpha/D;
		S(ii) = T_12/D;
	end
	T_out(2:N-1) = cumprod(hn.*S(1:N-2)).*S(2:N-1);
endif

T_out(1) = S(1);


endfunction	
