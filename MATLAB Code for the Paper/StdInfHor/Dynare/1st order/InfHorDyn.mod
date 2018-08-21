var w, Y, K, z, c, r;
predetermined_variables K;
varexo eps;
parameters calpha, cbeta, cdelta, cgamma, lbar, cmu, crho, csigma ;

calpha = .35;
cbeta = .96;
cdelta = .05;
cgamma = 2;
lbar = 1;
cmu = 0;
crho = .95^4;
csigma = 0.0130058759370760;


/*---------------------------------Model---------------------------------*/
model;
c = w*lbar + (1 + r - cdelta)*K - K(+1);
c^(-cgamma) = cbeta*((1 + r(+1) - cdelta)*(c(+1)^-cgamma));
Y = exp(z) * (K^calpha) * (lbar^(1 - calpha));
w = (1 - calpha) * exp(z) * (K/lbar)^calpha;
r = calpha * exp(z) * (lbar/K)^(1 - calpha);
z = crho * z(-1) + (1 - crho)*cmu + eps;
end; 

initval;
c	= .75;
r	=.12;
w	=1.32;
K	= 8.04;
z	= 0;
Y	= 2;
end;

steady;

/*-----------------------------Shock-------------------------------------*/

shocks;
var eps ;
stderr csigma ;
end;

/*---------------------------Simluation----------------------------------*/
stoch_simul(order = 1, irf = 150, nodisplay);
//Runs IRF func when we need irfs ^
//stoch_simul(order = 1, periods = 500);
//Runs stochastic simul when needed