%-------------------------------------------------------------------------%
% Set parameters for solutions to growth model with varying marginal tax
% rates
%
% This m-file calls the following m-file(s):
%     * tauchenhussey.m: Creates the discretized support for z and Markov
%          transition matrix for approximating integrals
%     * ksssolve.m: Solve for the certainty equivalent steady-state
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/varyingMtax/ParamSet') ;

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
% tau1      = scalar, asymptotic marginal tax rate as income goes to to
%             negative infinity
% tau2      = scalar, asymptotic marginal tax rate as income goes to to
%             positive infinity
% atau      = scalar, transition parameter in marginal tax rate function
% btau      = scalar, left/right shifter in marginal tax rate function
% parpi     = scalar, pi (to be used in dynare program)
%-------------------------------------------------------------------------%
beta       = 0.96 ;
alpha      = 0.35 ;
delta      = 0.05 ;
gamma      = 2.5 ;
L          = 1 ;

mu         = 0 ;
rho_qtr    = 0.95 ; %.95
rho        = rho_qtr^(4) ;
sigma_qtr  = 0.007 ; %.007
rhosigsum  = 0 ;
for sind = 1:4
    rhosigsum = rhosigsum + rho_qtr^(2*(sind-1)) ;
end
sigma = sqrt(rhosigsum)*sigma_qtr ;

Tirf  = 50 ;
Tsim  = 1000 ;

tau1 = 0.10 ; % 0.10
tau2 = 0.23 ; % 0.23
atau = 4 ; % 1, 2, 4, 8
btau = -12 ; % -3, -6, -12, -24

parpi = pi ;

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
% Solve for certainty equivalent steady-state at which the marginal tax
% rate is on the low and flat part of the arctan function
%-------------------------------------------------------------------------%
% kssinit = scalar, initial guess for steady state capital stock
% kss     = scalar, certainty-equivalent steady-state capital stock
% rss     = scalar, certainty-equivalent steady-state interest rate
% wss     = scalar, certainty-equivalent steady-state real wage
% css     = scalar, certainty-equivalent steady-state consumption
% tauss   = scalar, certainty-equivalent steady-state marginal tax rate
%-------------------------------------------------------------------------%
kssinit = 1 ;
params = [alpha,beta,delta,tau1,tau2,atau,btau] ;
options = optimset('Display','off','MaxFunEvals',100000,'MaxIter',1000,...
                    'TolFun',1e-15) ;
[kss,fval_kss] = fsolve(@ksssolve,kssinit,options,params) ;
display(kss) ;
display(fval_kss) ;

rss = alpha*(kss^(alpha-1)) ;
wss = (1-alpha)*(kss^alpha) ;
css = wss + (rss - delta)*kss ;
tauss = tau1 + ((1/pi)*atan(atau*(wss + (rss - delta)*kss) + btau) + 0.5)*...
        (tau2-tau1) ;
    
%-------------------------------------------------------------------------%
% Plot figure of tax function highlighting wss+(rss-delta)*kss and tauss
%-------------------------------------------------------------------------%
cmin = 10^(-6) ;
cmax = 7 ;
csize = 1000 ;
cvals = linspace(cmin,cmax,csize) ;
mtaxvals = tau1*ones(size(cvals)) + ...
           ((1/pi)*atan(atau*cvals + btau*ones(size(cvals))) + ...
           0.5*ones(size(cvals))).*(tau2 - tau1) ;

figure(1)
plot(cvals,mtaxvals,'LineWidth',3) ;
line([css css],[0.8*tau1 tauss]) ;
line([0 css],[tauss tauss]) ;
xlabel('Gross taxable income','FontSize', 12) ;
ylim([0.8*tau1,1.1*tau2])
ylabel('Marginal tax rate (\tau)','FontSize', 12) ;

%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save varyingMtax_par.mat ;
