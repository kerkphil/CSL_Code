%-------------------------------------------------------------------------%
% Set parameters for solutionsto standard infinite horizon problem for
% which the closed form solution is known
%
% This m-file calls the following m-file(s):
%     * tauchenhussey.m: Creates the discretized support for z and Markov
%          transition matrix for approximating integrals
%
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/ClosedForm') ;

starttime = clock ; % start timer

%-------------------------------------------------------------------------%
% Set parameters
%-------------------------------------------------------------------------%
% beta      = per-period discount factor
% alpha     = capital share of income, production function parameter
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
Tirf  = 40 ;
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
%epsvec4 = [0.001*sigma, zeros(1,Tirf-1)] ;
zirfvec1 = zeros(1, Tirf) ;
zirfvec2 = zeros(1, Tirf) ;
zsimvec3 = zeros(1,Tsim) ;
%zirfvec4 = zeros(1, Tirf) ;
for tind = 1:Tirf
    if tind == 1
        zirfvec1(tind) = mu + epsvec1(tind) ;
        zirfvec2(tind) = mu + epsvec2(tind) ;
        %zirfvec4(tind) = mu + epsvec4(tind) ;
    else
        zirfvec1(tind) = rho*zirfvec1(tind-1) + (1-rho)*mu + epsvec1(tind);
        zirfvec2(tind) = rho*zirfvec2(tind-1) + (1-rho)*mu + epsvec2(tind);
        %zirfvec4(tind) = rho*zirfvec4(tind-1) + (1-rho)*mu + epsvec4(tind);
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
% Solve for certainty equivalent steady-state
%-------------------------------------------------------------------------%
% rss = scalar, certainty-equivalent steady-state interest rate
% kss = scalar, certainty-equivalent steady-state capital stock
%-------------------------------------------------------------------------%
rss = 1/beta ;
kss = (alpha/rss)^(1/(1-alpha)) ;

%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save ClosedFormParams.mat ;