function e = Mod1dyn(x,params)
% calculates the deviation of the Euler equation from zero for input values
% of K(t-1), K(t), K(T+1), z(t) & z(t+1)
% note we want this to be zero in the steady state (not one) so that we can
% use fsolve to find SS values, if we wish.

beta = params(2) ;
delta = params(3) ;
gamma = params(4) ;

% read in values
Kplus  = x(1);
Know   = x(2);
Kminus = x(3);
zplus  = x(4);
znow   = x(5);

% get current period definitions
temp = Mod1defs(Kminus, znow, Know, params);
Ynow = temp(1);
wnow = temp(2);
rnow = temp(3);
cnow = temp(4);

% get next period definitions
temp  = Mod1defs(Know, zplus, Kplus, params);
Yplus = temp(1);
wplus = temp(2);
rplus = temp(3);
cplus = temp(4);

% calculate Euler equations
e = 1-beta*(1 + rplus - delta)*((cnow/cplus)^gamma) ;