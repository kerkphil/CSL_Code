%-------------------------------------------------------------------------%
% Solve simple RBC model by higher order Taylor series approximation using
% DYNARE
%-------------------------------------------------------------------------%
% This m-file calls the following m-file(s):
%-------------------------------------------------------------------------%

var y c k i r w z ;
varexo eps ;

parameters alpha beta delta lbar gamma mu sigma rho ;

alpha = 0.3 ;
beta  = 0.99 ;
delta = 0.025 ;
lbar  = 1 ;
gamma = 2 ;

mu    = 0 ;
sigma = 0.05 ;
rho   = 0.9 ;

model ;
c^(-gamma) = beta*(1 + r(+1) - delta)*(c(+1)^(-gamma)) ;
c + k = w*lbar + (1 + r - delta)*k(-1) ;
y = exp(z)*(k(-1)^alpha)*(lbar^(1-alpha)) ;
i = y - c ;
r = alpha*exp(z)*(lbar/k(-1))^(1-alpha) ;
w = (1-alpha)*exp(z)*((k(-1)/lbar)^alpha) ;
z = rho*z(-1) + (1 - rho)*mu + eps ;
end ;

initval ;
y = exp(mu)*(21.6333351205289^alpha)*(lbar^(1-alpha)) ;
c = ((1-alpha)*exp(mu)*(21.6333351205289/lbar)^alpha)*lbar +
    (1 + (alpha*exp(mu)*(lbar/21.6333351205289)^(1-alpha)) - delta)*
    21.6333351205289 ;
k = 21.6333351205289 ;
r = alpha*exp(mu)*(lbar/21.6333351205289)^(1-alpha) ;
w = (1-alpha)*exp(mu)*(21.6333351205289/lbar)^alpha ;
z = mu ;
end ;

steady(solve_algo=0) ; % make sure that I use the best method (0,1,2,3)

shocks ;
var eps;
stderr sigma;
end ;

stoch_simul(order=2, irf=200) ;


save TayApproxDyn.mat ;


