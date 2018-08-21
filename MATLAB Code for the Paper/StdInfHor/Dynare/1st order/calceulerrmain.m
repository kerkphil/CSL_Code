% calceulerrmain
%-------------------------------------------------------------------------%
% This file takes the series k and z to calc Euler Errors from the Value
% Function Iteration for a standard infinite horizon model.  It takes the
% equation and integrates over the conditional probabilities
%-------------------------------------------------------------------------%
% Requires that you pass in a value of ks and zs
%-------------------------------------------------------------------------%
% Euler error is eta_t = beta * E[(1 + r_t+1 - delta) * u'(c_t+1)] / u'(ct)

load('K_series.mat')
load('zshockdynirf1')
load('params')
load('kmatirf1ord.mat')
load('kmatsimseries1ord.mat')
load('kmat2irf1ord.mat')
load('simulationshocks.mat')


shock2eps = zeros(size(K_eps_dyn1)) ;
shock2eps(2) = csigma ;

% We are using the k chosen yesterday for today which we will call k_t and
% the z that was from yesterday z_t-1 and today's shock eps_t to calc k_t+1
z_eps = [0; z_eps] ;
params = [cbeta, cgamma, cdelta, mu, csigma, crho, calpha] ;
params2 = [ystar, A, B] ;

% Create space to fill euler errors
eedynirf1 = zeros(length(K_eps), 1) ;
K_eps = K_eps + ones(length(K_eps), 1) * kss ;

% Calculate Euler Errors for 1 st dev irf
for i = 1:length(K_eps) ;
    % We should be taking z_{t-1}, but we needed an initial z.  We set
    % z_{t-1} = 0 in period 0 and everything else starts in period 1
    eedynirf1(i) = calceulerrdyn1(K_eps(i), z_eps(i), shock2eps(i), params, params2) ;
end

RMSEirf = sqrt(sum(eedynirf1.^2)) ;
Maxeeirf = max(abs(eedynirf1)) ;
Mineeirf = min(abs(eedynirf1)) ;


% Create space to fill euler errors
eedynsersim1 = zeros(length(kmatsim), 1) ;

% Calculate Euler Errors for 1 st dev irf
for i = 1:length(kmatsim) ;
    % We should be taking z_{t-1}, but we needed an initial z.  We set
    % z_{t-1} = 0 in period 0 and everything else starts in period 1
    eedynsersim1(i) = calceulerrdyn1(kmatsim(i), zmatsim(i), simeps(i), params, params2) ;
end

RMSEsim = sqrt(sum(eedynsersim1.^2)) ;
Maxeesim = max(abs(eedynsersim1)) ;
Mineesim = min(abs(eedynsersim1)) ;


% Create space to fill euler errors
eedynnss = zeros(length(kmat2), 1) ;
zmat2 = [zmat2] ;

% Calculate Euler Errors for 1 st dev irf
for i = 1:length(kmat2) ;
    % We should be taking z_{t-1}, but we needed an initial z.  We set
    % z_{t-1} = 0 in period 0 and everything else starts in period 1
    eedynnss(i) = calceulerrdyn1(kmat2(i), zmat2(i), shock2eps2(i), params, params2) ;
end

RMSEnss = sqrt(sum(eedynnss.^2)) ;
Maxeenss = max(abs(eedynnss)) ;
Mineenss = min(abs(eedynnss)) ;

