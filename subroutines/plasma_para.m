% calculate the plasma angular frequency and relative permittivity at a frequency point.
% written by ZHU Bo at Nanjing University (bzhu@nju.edu.cn)

function [er_p wp] = plasma_para(ne, w)
% ne: electron density
% w: angular frequency

e0 = 8.854e-12;
me = 9.1e-31; % electron mass
e =  1.602e-19; % electron charge
wp = sqrt(ne/me/e0)*e;
er_p = 1 - (wp/w)^2;
