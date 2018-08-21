%-------------------------------------------------------------------------%
% Set parameters for solutions to growth model with three steady states
%
% This m-file calls the following m-file(s):
%     * tauchenhussey.m: Creates the discretized support for z and Markov
%          transition matrix for approximating integrals
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/ustab') ;

starttime = clock ; % start timer

%-------------------------------------------------------------------------%
% Set parameters
%-------------------------------------------------------------------------%
% beta      = per-period discount factor
% alpha     = capital share of income, production function parameter
% delta     = rate of depreciation
% gamma     = coefficient of relative risk aversion
% L         = individual endowment of inelastically supplied labor
% mu        = unconditional mean of AR(1) process for z
% rho_qtr   = quarterly persistence in AR(1) process for z: Hansen (1985)
% rho       = per-period persistence in AR(1) process for z
% sigma_qtr = quarterly standard deviation of shock in z process: Hansen
%             (1985)
% rhosigsum = coefficient one per-period standard deviation
% sind      = period index used in making sum
% sigma     = per-period standard deviation of shock in z process
% Tirf      = scalar, number of periods to simulate in IRFs
% Tsim      = scalar, number of periods to simulate in simulation test
%-------------------------------------------------------------------------%
beta       = 0.96 ;
alpha      = 0.35 ;
delta      = 0.05 ;
gamma      = 2.5 ;
L          = 1 ;
mu         = 0 ;
rho_qtr    = 0.95 ;
rho        = rho_qtr^(4) ;
sigma_qtr  = 0.007 ;
rhosigsum  = 0 ;
for sind = 1:4
    rhosigsum = rhosigsum + rho_qtr^(2*(sind-1)) ;
end
sigma = sqrt(rhosigsum)*sigma_qtr ;
Tirf  = 50 ;
Tsim  = 1000 ;

%-------------------------------------------------------------------------%
% Use Tauchen-Hussey method to generate Markov transition matrix
%-------------------------------------------------------------------------%
% znodes    = scalar, number of nodes in the discretized support of z used
%             in approximating expectations
% sigmaZ    = scalar, parameter used in the Tauchen method
% sigmaW    = scalar, parameter used in the Tauchen method
% baseSigma = scalar, parameter used in the Tauchen method
% zvec      = znodes x 1 vector, discretized support of z
% zTrans    = znodes x znodes matrix, Tauchen-Hussey Markov transition
%             matrix for z
%-------------------------------------------------------------------------%
znodes = 7 ;
sigmaZ = sigma/sqrt(1-rho^2) ;
sigmaW = 0.5 + rho/4 ;
baseSigma = sigmaW*sigma + (1-sigmaW)*sigmaZ ;
[zvec, zTrans] = tauchenhussey(znodes,mu,rho,sigma,baseSigma) ;

%-------------------------------------------------------------------------%
% Generate 3 different vectors of productivity levels corresponding to the
% two IRFs and the one simulation
%-------------------------------------------------------------------------%
% epsvec1  = 1 x Tirf vector of productivity shocks for IRF1, one standard
%            deviation shock. eps = [sigma,0,0,...0]
% epsvec2  = 1 x Tirf vector of productivity shocks for IRF2, two standard
%            deviation shock. eps = [2*sigma,0,0,...0]
% epsvec3  = 1 x Tsim vector of productivity shocks for the simulation
% zirfvec1 = 1 x Tirf vector of productivity levels for IRF1
% zirfvec2 = 1 x Tirf vector of productivity levels for IRF2
% zsimvec3 = 1 x Tsim vector of productivity levels for the simulation
%-------------------------------------------------------------------------%
epsvec1 = [sigma, zeros(1,Tirf-1)] ;
epsvec2 = [2*sigma, zeros(1,Tirf-1)] ;
epsvec3 = sigma*randn(1,Tsim) ;
zirfvec1 = zeros(1, Tirf) ;
zirfvec2 = zeros(1, Tirf) ;
zsimvec3 = zeros(1,Tsim) ;
for tind = 1:Tirf
    if tind == 1
        zirfvec1(tind) = mu + epsvec1(tind) ;
        zirfvec2(tind) = mu + epsvec2(tind) ;
    else
        zirfvec1(tind) = rho*zirfvec1(tind-1) + (1-rho)*mu + epsvec1(tind);
        zirfvec2(tind) = rho*zirfvec2(tind-1) + (1-rho)*mu + epsvec2(tind);
    end
end
for tind = 1:Tsim
    if tind == 1
        zsimvec3(tind) = mu + epsvec3(tind) ;
    else
        zsimvec3(tind) = rho*zsimvec3(tind-1) + (1-rho)*mu + epsvec3(tind);
    end
end

%-------------------------------------------------------------------------%
% Solve for certainty-equivalent unstable steady-state
%-------------------------------------------------------------------------%
% rss   = scalar, certainty-equivalent steady-state interest rate
% kss   = scalar, certainty-equivalent steady-state capital stock
% wss   = scalar, certainty-equivalent steady-state real wage
% css   = scalar, certainty-equivalent steady-state consumption
% d     = scalar, fixed lump sum transfer from government
% tauss = scalar, certainty-equivalent steady-state marginal tax rate
% Bss   = scalar, certainty-equivalent steady-state government bond balance
%-------------------------------------------------------------------------%
rss = 1/beta + delta - 1 ;
kss = (alpha/rss)^(1/(1-alpha)) ;
wss = (1-alpha)*(kss^alpha) ;
d = 0.05*kss ;
tauss = d/wss ;
css = (1 - tauss)*wss + (rss - delta)*kss + d ;
Bss = 0 ;

%-------------------------------------------------------------------------%
% Solve for certainty-equivalent stable low steady state
%-------------------------------------------------------------------------%
% Bmin   = scalar, minimum cutoff value of government debt, below which the
%          marginal tax rate is 0.99
% Blow   = scalar, low cutoff value of government debt. The marginal tax
%          rate between Bmin and Blow is a linear function between tauss
%          and 0.99
% BLss   = scalar, low government debt certainty-equivalent steady-state
%          government debt level
% tauLss = scalar, low government debt certainty-equivalent steady-state
%          marginal tax rate
% kLss   = scalar, low government debt certainty-equivalent steady-state
%          individual capital holdings
% cLss   = scalar, low government debt certainty-equivalent steady-state
%          individual consumption
%-------------------------------------------------------------------------%
Bmin = -kss ;
Blow = -0.5*kss ;
BLss = (((0.99 - tauss)/(Bmin - Blow))*Blow + (d/wss) - tauss)/...
       ((0.99 - tauss)/(Bmin - Blow) + (rss - delta)/wss) ;
tauLss = ((0.99 - tauss)/(Bmin - Blow))*BLss + tauss - ...
         ((0.99 - tauss)/(Bmin - Blow))*Blow ;
kLss = (alpha/rss)^(1/(1 - alpha)) - BLss ;
cLss = (1 - tauLss)*wss + (rss - delta)*kLss + d ;

%-------------------------------------------------------------------------%
% Solve for certainty-equivalent stable high steady state
%-------------------------------------------------------------------------%
% Bmax   = scalar, maximum cutoff value of government savings, above which
%          the marginal tax rate is zero
% Bhigh  = scalar, high cutoff value of government savings. The marginal
%          tax rate between Bhigh and Bmax is a linear function between
%          tauss and 0
% BHss   = scalar, high government savings certainty-equivalent steady-
%          state government debt level
% tauHss = scalar, high government savings certainty-equivalent steady-
%          state marginal tax rate
% kHss   = scalar, high government savings certainty-equivalent steady-
%          state individual capital holdings
% cHss   = scalar, high government savings certainty-equivalent steady-
%          state individual consumption
%-------------------------------------------------------------------------%
Bmax = kss ;
Bhigh = 0.5*kss ;
BHss = ((tauss/(Bhigh - Bmax))*Bhigh + (d/wss) - tauss)/...
       (tauss/(Bhigh - Bmax) + (rss - delta)/wss) ;
tauHss = (tauss/(Bhigh - Bmax))*BHss + tauss - ...
         Bhigh*(tauss/(Bhigh - Bmax)) ;
kHss = (alpha/rss)^(1/(1 - alpha)) - BHss ;
cHss = (1 - tauHss)*wss + (rss - delta)*kHss + d ;

%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save unstab_par.mat ;
