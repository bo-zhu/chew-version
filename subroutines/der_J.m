% calculate der_\hat{J}(x), where \hat{J}(x) = x*sbesselj(x),
% and der_\hat{J}(x) is the 1st order derivative of \hat{J}(x).
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [Y] = der_J(ALPHA,X,OPT)
% ALPHA: degree of the functions.
% OPT: if available, return the normalized value.

Y = sbesselj(ALPHA,X,OPT)+X*der_sbesselj(ALPHA,X,OPT) ;

endfunction
