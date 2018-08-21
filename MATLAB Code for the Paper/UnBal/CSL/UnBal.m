% Unbalanced Growth Model
clear

% set model parameters
alf = .35;
bet = .98;
sig = .013;
rho = .95;
gam = 2.5;
del = .025;
g   = .01;
chi = 5;
theta = 1.5;

% set starting values
X0 = [10; .5; 1];

% set up parameter vector to pass to DSGE function file
param = [alf bet sig rho gam del g chi theta];
% set numerical parameters
nx = 3;
ny = 0;
nz = 1;
nobs = 1000;
randomerr = 1;
logX = 0;

% generate discret support for epsilon to be used in Euler error calcs
ne = 100;  %number of elements in support
Eps = zeros(ne,1);
Cum = -.5/ne;
Phi = ones(ne,1)/ne;
for e = 1:ne
    Cum = Cum + Phi(e);
    Eps(e) = norminv(Cum,0,sig);
end

% generate a history of Z's
Z = zeros(nobs+2,nz);
if randomerr
    eps = randn(nobs+2,1)*sig;
else
    load simulationshocks.mat
    eps = simeps;
    [nobs,~] = size(simeps)-2;
end
for t=1:nobs+1
    Z(t+1,:) = Z(t,:)*rho + eps(t+1,:);
end

%  current state linarization
tic;
% [XCSL, ~] = LinApp_CSL(@UnBal_dyn,param,X0',Z,rho,logX);
[XCSL, temp, EulerErr] = LinApp_CSL_Euler(@UnBal_dyn,param,X0',Z,...
    rho,logX,Eps,Phi);
k  = XCSL(:,1);
ell  = XCSL(:,2);
Y = zeros(nobs+2,1);
r  = zeros(nobs+2,1);
w  = zeros(nobs+2,1);
c = zeros(nobs+2,1);
i  = zeros(nobs+2,1);

for t=2:nobs+1
    [Y(t+1), r(t+1), w(t+1), c(t+1), i(t+1)] = ...
     UnBal_defs(k(t), ell(t), t, Z(t), k(t+1), ell(t+1), t+1, param);
end
toc

AvgEE = mean(EulerErr(2:nobs+1,:))
MaxAEE = max(abs(EulerErr(2:nobs+1,:)))
RMSEE = sqrt(mean(EulerErr(2:nobs+1,:).^2))

%  plot simulation
UnBal_Plot1(k(2:nobs+1,:),ell(2:nobs+1,:))
