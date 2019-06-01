% First derivative of spherical Hankel functions of the first (K=1) or second (K=2) kind
% Written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [DER_H] = der_sbesselh(ALPHA, K, X, OPT)
% ALPHA: degree of the function.
% OPT: if available, return the normalized value.

DER_H = ALPHA./X.*sbesselh(ALPHA,K,X,OPT) - sbesselh(ALPHA+1,K,X,OPT);

endfunction
