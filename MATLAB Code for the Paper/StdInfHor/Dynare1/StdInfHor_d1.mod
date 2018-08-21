var w, k, z, c, r;
/*predetermined_variables K;*/
varexo eps;
parameters calpha, cbeta, cdelta, cgamma cmu, crho, csigma ;

load '../ParamSet/StdInfHorParams.mat' alpha beta delta gamma mu rho sigma ;

calpha = alpha ;
cbeta = beta ;
cdelta = delta ;
cgamma = gamma ;
cmu = mu ;
crho = rho ;
csigma = sigma ;

clear alpha beta delta gamma mu rho sigma T ;

/*---------------------------------Model---------------------------------*/
model ;
    c ^ (-cgamma) = cbeta * (1 + r(+1) - cdelta) * c(+1) ^ (-cgamma) ;
    c = w + (1 + r - cdelta) * k(-1) - k ;
    w = (1-calpha) * exp(z) * k(-1) ^ calpha ;
    r = calpha * exp(z) * (k(-1) ^ (calpha - 1)) ;
    z = crho * z(-1) + (1 - crho) * cmu + eps ;
end ; 

initval ;
    c = 1.665;
    r = 0.0917;
    w = 1.337;
    k = 7.855;
    z = cmu ;
end ;

steady ;

/*-----------------------------Shock-------------------------------------*/

shocks;
var eps ;
stderr csigma ;
end;

/*---------------------------Simluation----------------------------------*/
stoch_simul(order = 1, irf = 50 ) ;