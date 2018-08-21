var w, k, z, c, r;
/*predetermined_variables K;*/
varexo eps;
parameters calpha, cbeta, cmu, crho, csigma ;

load '../ParamSet/ClosedFormParams.mat' alpha beta mu rho sigma ;

calpha = alpha ;
cbeta = beta ;
cmu = mu ;
crho = rho ;
csigma = sigma ;

clear alpha beta mu rho sigma T ;

/*---------------------------------Model---------------------------------*/
model ;
    (1 / c) = cbeta * r(+1) / c(+1) ;
    c = w + r*k(-1) - k ;
    w = (1-calpha) * exp(z) * k(-1) ^ calpha ;
    r = calpha * exp(z) * (k(-1) ^ (calpha - 1)) ;
    z = crho * z(-1) + (1 - crho) * cmu + eps ;
end ; 

initval ;
    c = 0.37;
    r = 1.04;
    w = 0.36;
    k = 0.19;
    z = cmu ;
end ;

steady ;

/*-----------------------------Shock-------------------------------------*/

shocks;
var eps ;
stderr csigma ;
end;

/*---------------------------Simluation----------------------------------*/
stoch_simul(order = 2, irf = 40 ) ;