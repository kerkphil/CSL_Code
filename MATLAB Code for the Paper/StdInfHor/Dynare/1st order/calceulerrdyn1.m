function eulerr = calceulerrdyn1(kt, ztm1, epst, params, params2)
%-------------------------------------------------------------------------%
% This file takes the series k and z to calc Euler Errors from the Value
% Function Iteration for a standard infinite horizon model.  It takes the
% equation and integrates over the conditional probabilities
%-------------------------------------------------------------------------%
% Requires that you pass in a value of ks and zs
%-------------------------------------------------------------------------%
% Euler error is eta_t = beta * E[(1 + r_t+1 - delta) * u'(c_t+1)] / u'(ct)



% Import parameters necessary for calculations
beta  = params(1) ;
gamma = params(2) ;
delta = params(3) ;
mu    = params(4) ;
sigma = params(5) ;
rho   = params(6) ;
calpha = params(7) ;
ystar = params2(:, 1) ;
A     = params2(:, 2:3) ;
B     = params2(:, 4);

kss = ystar(3);

% Set up a group of functions that will calculate the value of individual
% variables given k and z
calczp = @(ztm, eps_shock) (rho*ztm + eps_shock) ; 
calcr  = @(ktt, ztt) (calpha .* exp(ztt) .*  ((ktt) ^ (calpha - 1))) ;
calcw  = @(ktt, ztt) ((1 - calpha) .* exp(ztt) .*  (ktt .^ calpha)) ;
calckp = @(kztempp, eps_tt) (ones(1, length(eps_tt)) .* (kss + A(3,:)*kztempp)...
    + B(3).*eps_tt) ;

% These are defs that will be used to calc other definitions
k_eps = kt - kss ;
ztm_eps = ztm1 - 0 ;
kztemp = [k_eps, ztm_eps]' ;

% Solve for necessary pieces of the euler equation
zt      = calczp(ztm1, epst) ;
L       = 1 ;
wt      = (1 - calpha) * exp(zt) * (kt/L) ^ calpha ;
rt      = calpha .* exp(zt) .* (L / kt) ^ (1 - calpha) ;
kprime  = calckp(kztemp, epst) ;
ct      = wt + (1 + rt - delta)*kt - kprime ;
muzt    = rho*zt + (1-rho)*mu ;

% These are defs that will be used to calc other definitions
kp_eps = kprime - kss ;
zt_eps = zt - 0 ;
kzptemp = [kp_eps, zt_eps]' ;


% Numer is the expected value part of euler error.  We will integrate over
% Numer to get expected value
% I think we want to integrate over the eps shock b/c of how Dynare
% calculates the policy functions.  Double check.
% Also, I have a problem with the policy function.  Why would our decision
% of K be dependent on today's shock?  That seems stupid.
Numer  = @(eps_t)pdf('norm',eps_t, 0, sigma) .* ...
    (((1 - delta) * ones(1, length(eps_t)) + calcr(kprime, calczp(zt, eps_t))) .* ...
    (calcw(kprime, calczp(zt, eps_t)) + ((1 - delta) * ones(1, length(eps_t)) + ...
    calcr(kprime, calczp(zt, eps_t))) .* kprime - ...
    calckp(kzptemp, eps_t)) .^ (-gamma)) ;

Denom  = ct ^ - gamma ;

Expect =  quadgk(Numer, -Inf, Inf) ;

eulerr = ((beta * Expect) / Denom) - 1 ;

return
