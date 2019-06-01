% First derivative of spherical Bessel functions of the first kind
% Written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [DER_J] = der_sbesselj(ALPHA, X, OPT)
% ALPHA: degree of the function.
% OPT: if available, return the normalized value.

DER_J = ALPHA./X.*sbesselj(ALPHA,X,OPT) - sbesselj(ALPHA+1,X,OPT);

endfunction
