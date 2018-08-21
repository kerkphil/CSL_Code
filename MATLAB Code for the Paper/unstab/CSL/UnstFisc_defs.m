function [y, i, c, r, w, tau, K] = UnstFisc_defs(k,B,z,kp,Bp,param)

alf = param(1);
del = param(6);
ell = param(7);
taufix = param(8);
d = param(9);
Bmax = param(10);
Bupp = param(11);
Blow = param(12);
Bmin = param(13);


K = k + B;
%Brat = B/K;
y = exp(z)*K^alf*ell^(1-alf);
i = kp - (1-del)*k;
r = alf*y/K;
w = (1-alf)*y/ell;
if B>Bmax
    tau = 0;
elseif B>Bupp
    tau = taufix*(Bmax-B)/(Bmax-Bupp);
elseif B>Blow
    tau = taufix;
elseif B>Bmin
    tau = taufix + (1-taufix)*(Blow-B)/(Blow-Bmin);
else
    tau = 1;
end
c = (1-tau)*w*ell+(1+r-del)*k - kp - d;
end


