% Spherical Hankel functions of the first (K=1) or second (K=2) kind.
% Written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [H IERR] = sbesselh(ALPHA, K, X, OPT)
% ALPHA: degree of the function.
% OPT: if availabel, return the normalized sbesselh.

if OPT=='norm'
	[H IERR] = besselh(ALPHA+0.5, K, X, OPT);
else 
	[H IERR] = besselh(ALPHA+0.5, K, X);
end

H = sqrt(pi./(2*X)).*H;

endfunction
