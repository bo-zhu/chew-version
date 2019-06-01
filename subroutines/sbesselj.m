% Spherical Bessel functions of the first kind
% Written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [J IERR] = sbesselj(ALPHA, X, OPT)
% ALPHA: degree of the function.
% OPT: if available, return the normalized value.

if OPT=='norm'
	[J IERR] = besselj(ALPHA+0.5, X, OPT);
else
	[J IERR] = besselj(ALPHA+0.5, X);
end

J = sqrt(pi./(2*X)).*J;


endfunction
