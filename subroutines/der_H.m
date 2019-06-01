% calculate der_\hat{H}(x), where \hat{H}(x) = x*sbesselh(x), 
% and der_\hat{H}(x) is the 1st order derivative of \hat{H}(x).
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [Y] = der_H(ALPHA,X,OPT)
% ALPHA: degree of the functions.
% OPT: if available, return the normalized value.

Y = sbesselh(ALPHA,1,X,OPT)+X*der_sbesselh(ALPHA,1,X,OPT) ;

endfunction
