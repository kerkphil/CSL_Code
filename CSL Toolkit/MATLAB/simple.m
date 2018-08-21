% Hansen's model without labor/leisure decision
clear

%set model parameters
A = 1;
theta = .33;
del = .025;
bet = .995;
hbar = .3;
gam = 1;
rho = .9;
sig = .02;
% set up parameter vector to pass to DSGE function file
param = [A theta del bet hbar gam rho sig];

%set numerical parameters
nx = 1;
ny = 0;
nz = 1;
options = optimset('Display','iter');
nobs = 250;
logX = true;
EE = 1;

% generate discret support for epsilon to be used in Euler error calcs
ne = 100;  %number of elements in support
Eps = zeros(ne,1);
Cum = -.5/ne;
Phi = ones(ne,1)/ne;
for e = 1:ne
    Cum = Cum + Phi(e);
    Eps(e) = norminv(Cum,0,sig);
end

Zbar = 0;
% find SS numerically
guessXY = 1;
Xbar = LinApp_FindSS(@simple_dyn,param,guessXY,Zbar,nx,ny);
theta0 = [Xbar; Xbar; Xbar; Zbar; Zbar];

NN = 0;
%find derivatives and coefficients numerically

[AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT] = ...
    LinApp_Deriv(@simple_dyn,param,theta0,nx,ny,nz,logX);

[PP, QQ, UU, RR, SS, VV] = ...
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar);

%find quadratic approximation's steady state\
%[XQtil1, XQtil2] = QuadRoots(.5*HXX,HX-1,H0+sig^2*Hvv/2);

%generate a history of Z's
Z = zeros(nobs,nz);
% uncomment for simulation
eps = sig*randn(nobs,nz);
% uncomment for IRF
% eps = zeros(nobs,nz);
% eps(3,1) = sig;
for t=1:nobs-1
    Z(t+1,:) = Z(t,:)*NN + eps(t+1,:);
end

% set starting values and simulate
XYbar = Xbar;
X0 = Xbar;
tic;
[XSSL, ~] = LinApp_SSL(X0',Z,XYbar',NN,PP,QQ,UU,logX,EE,Eps,Phi,...
                       @simple_dyn,param);
toc

tic;
[XCSL, ~] = LinApp_CSL(@simple_dyn,param,X0',Z,NN,logX,EE,Eps,Phi);
toc

% plot results
plotdata = [XSSL XCSL];
figure;
plot(plotdata)
