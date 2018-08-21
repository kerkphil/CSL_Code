function [Y, r, w, c, i] = ...
    UnBal_defs(km, ellm, taum, z, k, ell, tau, param)

alf = param(1);
bet = param(2);
sig = param(3);
rho = param(4);
gam = param(5);
del = param(6);
g   = param(7); 
chi = param(8);
theta = param(9);

Y = km^alf*((1-ell)*exp(g*tau+z))^(1-alf);
r = alf*Y/km;
w = (1-alf)*Y/(1-ell);
c = w*(1-ell) + (1+r-del)*km - k;
i = k - (1-del)*km;