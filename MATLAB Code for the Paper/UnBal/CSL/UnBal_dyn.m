function out = UnBal_dyn(in,param)

kp    = in(1);
ellp  = in(2);
taup  = in(3);
k     = in(4);
ell   = in(5);
tau   = in(6);
km    = in(7);
ellm  = in(8);
taum  = in(9);
zp    = in(10);
z     = in(11);

alf = param(1);
bet = param(2);
sig = param(3);
rho = param(4);
gam = param(5);
del = param(6);
g   = param(7); 
chi = param(8);
theta = param(9);

[Y, r, w, c, i] = ...
    UnBal_defs(km, ellm, taum, z, k, ell, tau, param);
[Yp, rp, wp, cp, ip] = ...
    UnBal_defs(k, ell, tau, zp, kp, ellp, taup, param);

out1 = (chi*(1-ell)^theta)/(c^(-gam)) - 1;
out2 = bet*(cp/c)^(-gam)*(1+rp-del) - 1;
out3 = tau - taum - 1;
out = [out1; out2; out3];