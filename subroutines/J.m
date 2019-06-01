% calculate \hat{J}(x), where \hat{J}(x) = x*sbesselj(x).
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [Y] = J(ALPHA,X,OPT)
% ALPHA: degree of the functions.
% OPT: if available, return the normalized value.

Y = X*sbesselj(ALPHA,X,OPT) ;

endfunction
