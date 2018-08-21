function y = CSLStdInfHordefs(K, z, Kplus, params)
% calculates the values of Y(t), r(t), w(t), C(t), tau(t) and tax(t) for
% the input values of K(t), K(t+1) & z(t)

alpha = params(1) ;
delta = params(3) ;

Y = exp(z)*(K^alpha) ;
r = alpha*exp(z)*(K^(alpha - 1)) ;
w = (1 - alpha)*exp(z)*(K^alpha) ;
c = w  + (1 + r - delta)*K - Kplus ;

y = [Y; w; r; c];