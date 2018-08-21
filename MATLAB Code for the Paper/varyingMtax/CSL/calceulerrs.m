function eulerr = calceulerrs(kt, zt, UUt, params)

% This function calculates the value of the Euler errors for StdInfHor model
% Euler error is defined as 
% eta_it = beta * E[(1 + r_t+1 -delta)*u'(c_i,t+1)]/u'(c_i,t)

% Returns:
% eulerr
% Import parameters necessary for calculations

beta  = params(1) ;
gamma = params(2) ;
delta = params(3) ;
mu    = params(4) ;
sigma = params(5) ;
rho   = params(6) ;
calpha = params(7) ;
kss   = params(8) ;

% Set up a group of functions that will calculate the value of individual
% variables given k and z
calcr  = @(ktt, ztt) (calpha .* exp(ztt) .*  ((ktt) ^ (calpha - 1))) ;
calcw  = @(ktt, ztt) ((1 - calpha) .* exp(ztt) .*  (ktt .^ calpha)) ;
calckp = @(UUt, Xtm, Zt) Xtm.*exp(UUt') ;

% Solve for necessary pieces of the euler equation
L       = 1 ;
wt      = (1 - calpha) * exp(zt) * (kt/L) ^ calpha ;
rt      = calpha .* exp(zt) .* (L / kt) ^ (1 - calpha) ;
kprime  = calckp(UUt, kt, zt) ;
ct      = wt + (1 + rt - delta)*kt - kprime ;
muzt    = rho*zt + (1-rho)*mu ;

% Numer is the expected value part of euler error.  We will integrate over
% Numer to get expected value
Numer  = @(zz)pdf('norm',zz, mu, sigma) .* ...
    (((1 - delta) * ones(1, length(zz)) + calcr(kprime, zz)) .* ...
    (calcw(kprime, zz) + ((1 - delta) * ones(1, length(zz)) + ...
    calcr(kprime, zz)) .* kprime - ...
    calckp(UUt, kprime, zz)) .^ (-gamma)) ;

Denom  = ct ^ - gamma ;

Expect =  quadgk(Numer, -Inf, Inf) ;

eulerr = ((beta * Expect) / Denom) - 1 ;

return
