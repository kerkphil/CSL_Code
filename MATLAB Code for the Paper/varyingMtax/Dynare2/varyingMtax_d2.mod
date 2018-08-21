var w, k, z, c, r, tau, d ;
/*predetermined_variables K;*/
varexo eps;
parameters calpha, cbeta, cdelta, cgamma cmu, crho, csigma, ctau1, ctau2, ca, cb cpi ;

load '../ParamSet/varyingMtax_par.mat' alpha beta delta gamma mu rho sigma tau1 tau2 atau btau css rss wss tauss parpi ;

calpha = alpha ;
cbeta = beta ;
cdelta = delta ;
cgamma = gamma ;
cmu = mu ;
crho = rho ;
csigma = sigma ;
ctau1 = tau1 ;
ctau2 = tau2 ;
ca = atau ;
cb = btau ;
cpi = parpi ;

clear alpha beta delta gamma mu rho sigma T tau1 tau2 atau btau parpi ;

/*---------------------------------Model---------------------------------*/
model ;
    c ^ (-cgamma) = cbeta * (c(+1) ^ (-cgamma)) * (1 + (1 - tau(+1)) * (r(+1) - cdelta) - ((ctau2 - ctau1) / cpi) * ((ca * (r(+1) - cdelta)) / (1 + (ca * (w(+1) + (r(+1) - cdelta) * k) + cb)^2)) * (w(+1) + (r(+1) - cdelta) * k))  ;
    c = w + (1 + r - cdelta) * k(-1) - tau * (w + (r - cdelta) * k(-1)) - k + d ;
    w = (1-calpha) * exp(z) * k(-1) ^ calpha ;
    r = calpha * exp(z) * (k(-1) ^ (calpha - 1)) ;
    tau = ctau1 + ((1/cpi) * atan(ca * (w + (r - cdelta) * k(-1)) + cb) + 0.5) * (ctau2 - ctau1) ;
    d = tau * (w + (r - cdelta) * k(-1)) ;
    z = crho * z(-1) + (1 - crho) * cmu + eps ;
end ; 

initval ;
    c = 1.2961 ;
    r = 0.0971 ;
    w = 1.2961 ;
    k = 7.1838 ;
    z = cmu ;
    tau = 0.1075 ;
    d = 0.1757 ;
end ;

steady ;

/*-----------------------------Shock-------------------------------------*/

shocks;
var eps ;
stderr csigma ;
end;

/*---------------------------Simluation----------------------------------*/
stoch_simul(order = 2, irf = 50 ) ;