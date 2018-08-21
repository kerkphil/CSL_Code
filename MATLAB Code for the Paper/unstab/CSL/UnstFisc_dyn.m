function out = UnstFisc_dyn(in,param)

bpp =  in(1);
Bpp =  in(2);
bp =   in(3);
Bp =   in(4);
b =    in(5);
B =    in(6);
zp =   in(7);
z =    in(8);

bet = param(2);
gam = param(5);
del = param(6);
ell = param(7);
d = param(9);

[~, ~, c, r, w, tau, ~] = UnstFisc_defs(b,B,z,bp,Bp,param);
[~, ~, cp, rp, ~, ~, ~] = UnstFisc_defs(bp,Bp,zp,bpp,Bpp,param);

out1 = bet*(c/cp)^gam*(1+rp-del) - 1;
out2 = Bp - tau*w*ell - (1+r-del)*B + d;
out = [out1; out2];