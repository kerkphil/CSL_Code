function e = CSLvaryingMtaxdyn(x,params)
% calculates the deviation of the Euler equation from zero for input values
% of K(t-1), K(t), K(T+1), z(t) & z(t+1)
% note we want this to be zero in the steady state (not one) so that we can
% use fsolve to find SS values, if we wish.

beta = params(2) ;
delta = params(3) ;
gamma = params(4) ;
tau1 = params(12);
tau2 = params(13);
a = params(14);
b = params(15);

% read in values
kplus  = x(1);
know   = x(2);
kminus = x(3);
zplus  = x(4);
znow   = x(5);

% get current period definitions
temp = CSLvaryingMtaxdefs(kminus, znow, know, params);
cnow = temp(6);

% get next period definitions
temp  = CSLvaryingMtaxdefs(know, zplus, kplus, params);
wplus = temp(2);
rplus = temp(3);
tauplus = temp(4);
cplus = temp(6);

% calculate Euler equations
deriv = ((tau2-tau1)/pi()) * ((a * (rplus - delta)) / ...
    (1 + (a*(wplus+(rplus-delta)*kplus)+b)^2)); 
e = cnow^(-gamma) - beta * ...
    (cplus^(-gamma)) * (1+(1-tauplus)*(rplus-delta) - ...
     deriv*(wplus+(rplus-delta)*kplus));