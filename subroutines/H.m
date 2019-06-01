% calculate \hat{H}(x), where \hat{H}(x) = x*sbesselh(x).
% written by ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )  

function [Y] = H(ALPHA,X,OPT)
% ALPHA: degree of the functions.
% OPT: if available, return the normalized value.

Y = X*sbesselh(ALPHA,1,X,OPT);

endfunction
