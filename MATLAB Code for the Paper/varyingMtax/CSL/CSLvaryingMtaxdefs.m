function out = CSLvaryingMtaxdefs(k, z, kp, params)
% calculates the values of y(t), r(t), w(t), c(t), tau(t) and d(t) for
% the input values of k(t), k(t+1) & z(t)

alpha = params(1) ;
delta = params(3) ;
tau1 = params(12);
tau2 = params(13);
a = params(14);
b = params(15);

y = exp(z)*k^alpha ;
w = (1-alpha)*y ;
r = alpha*y/k ;
tau = tau1 + ((1/pi()) * atan(a*(w+(r-delta)*k)+b)+0.5) * (tau2-tau1) ;
d = tau * (w + (r - delta) * k) ;
c = w + (1 + r - delta) * k - tau * (w + (r - delta) * k) - kp + d ;

out = [y; w; r; tau; d; c];