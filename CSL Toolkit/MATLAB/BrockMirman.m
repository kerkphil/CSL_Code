% Borck & Mirman model
clear

%set model parameters
alf = .35;
bet = .98;
sig = .02;
rho = .95;
% set up parameter vector to pass to DSGE function file
param = [alf bet sig rho];

%set numerical parameters
nx = 1;
ny = 0;
nz = 1;
options = optimset('Display','iter');
nobs = 250;
logX = 0;
DO_QZ = 0;
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
XYbar = LinApp_FindSS(@BrockMirman_dyn,param,.1,Zbar,nx,ny);
Xbar = XYbar(1:nx);
Ybar = XYbar(nx+1:nx+ny);
theta0 = [Xbar; Xbar; Xbar; Zbar; Zbar];

NN = rho;
%find derivatives and coefficients numerically

[AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT] = ...
    LinApp_Deriv(@BrockMirman_dyn,param,theta0,nx,ny,nz,logX);

[PP, QQ, UU, RR, SS, VV] = ...
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar)

%generate a history of Z's
Z = zeros(nobs,nz);
eps = sig*randn(nobs,nz);
for t=1:nobs-1
    Z(t+1,:) = Z(t,:)*NN + eps(t+1,:);
end

% set starting values and simulate
XYbar = Xbar;
X0 = Xbar;

%  steady state linarization
tic;
[XSSL, ~] = LinApp_SSL(X0',Z,XYbar',NN,PP,QQ,UU,logX,EE,Eps,Phi,...
                       @BrockMirman_dyn,param);
toc

%  current state linarization
tic;
[XCSL, ~] = LinApp_CSL(@BrockMirman_dyn,param,X0',Z,NN,logX,EE,Eps,Phi);
toc

%  exact solution
Xexact = zeros(nobs,nx);
Xexact(1,:) = X0;
for t=1:nobs-1
    Xexact(t+1,:) = alf*bet*exp(Z(t+1,:))*Xexact(t,:)^alf;
end

plotdata = [XSSL XCSL Xexact];
ratiodata = [log(XSSL./Xexact) log(XCSL./Xexact)];

figure;
plot(plotdata) 
figure;
plot(ratiodata)
