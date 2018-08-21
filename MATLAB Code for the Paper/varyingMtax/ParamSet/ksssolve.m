function diff = ksssolve(kval,params)

alpha = params(1) ;
beta = params(2) ;
delta = params(3) ;
tau1 = params(4) ;
tau2 = params(5) ;
a = params(6) ;
b = params(7) ;

%--------------------------------------------------------------------------
% Calculate the error (difference) in the steady-state certainty-equivalent
% equilibrium household first order condition
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
r = alpha*(kval^(alpha-1)) ;
w = (1 - alpha)*(kval^alpha) ;
x = w + (r - delta)*kval ;
tau = tau1 + ((1/pi)*atan(a*x + b) + 0.5)*(tau2 - tau1) ;
dtaudk = ((tau2 - tau1)/pi)*((a*alpha*(kval^(alpha - 1)) - a*delta)/...
         (1 + (a*x + b)^2)) ;

diff = 1 + (1 - tau)*(r - delta) - dtaudk*(w + (r - delta)*kval) - 1/beta ;
