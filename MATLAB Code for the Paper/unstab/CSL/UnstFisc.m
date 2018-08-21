% Unstable Fiscal Policy Model
clear

% set model parameters
alf = .35;
bet = .98;
sig = 0.013;
rho = .95;
gam = 2.0;
del = .025;
ell = 1;
taufix = .20;
d = .5;
Bmax = 50;
Bupp = 10;
Blow = -10;
Bmin = -50;

% set up parameter vector to pass to DSGE function file
param = [alf bet sig rho gam del ell taufix d Bmax Bupp Blow Bmin];
% set numerical parameters
nx = 2;
ny = 0;
nz = 1;
options = optimset('Display','iter');
nobs = 500;
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

% find SS numerically
Zbar = 0;
guess = [1; .01];
XYbar = LinApp_FindSS(@UnstFisc_dyn,param,guess,Zbar,nx,ny);
Xbar = XYbar(1:nx)
Ybar = XYbar(nx+1:nx+ny);
theta0 = [Xbar; Xbar; Xbar; Zbar; Zbar];
check = UnstFisc_dyn(theta0,param)
kbar = Xbar(1);
Bbar = Xbar(2);
% find other SS values
[ybar, ibar, cbar, rbar, wbar, taubar, Kbar] =...
    UnstFisc_defs(kbar,Bbar,Zbar,kbar,Bbar,param);


% find derivatives and coefficients numerically
NN = rho;
[AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT] = ...
    LinApp_Deriv(@UnstFisc_dyn,param,theta0,nx,ny,nz,logX);
[PP, QQ, UU, RR, SS, VV] = ...
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,Zbar);


% generate a history of Z's
Z = zeros(nobs+2,nz);
if randomerr
    eps = randn(nobs+2,1)*sig;
else
    load simulationshocks.mat
    eps = simeps;
    [nobs,~] = size(simeps);
    nobs = nobs-2;
end
for t=1:nobs+1
    Z(t+1,:) = Z(t,:)*NN + eps(t+1,:);
end


% set starting values and simulate
X0 = Xbar;
X0(2) = 1*Xbar(2);

%  steady state linarization
tic;
[XSSL, temp, EulerErrSSL] = LinApp_SSL_Euler(X0',Z,Xbar',logX,...
    PP,QQ,UU,NN,Eps,Phi,param,@UnstFisc_dyn);
kSSL = XSSL(:,1);
BSSL = XSSL(:,2);
ySSL = zeros(nobs+2,1);
iSSL = zeros(nobs+2,1);
cSSL = zeros(nobs+2,1);
rSSL = zeros(nobs+2,1);
wSSL = zeros(nobs+2,1);
tauSSL = zeros(nobs+2,1);
KSSL = zeros(nobs+2,1);
for t=1:nobs+1
    [ySSL(t+1),iSSL(t+1),cSSL(t+1),rSSL(t+1),wSSL(t+1),tauSSL(t+1)...
        ,KSSL(t+1)] = UnstFisc_defs(kSSL(t),BSSL(t),Z(t+1),kSSL(t+1),...
        BSSL(t+1),param);
end

AvgEESSL = mean(EulerErrSSL(2:nobs+1,:))
MaxAEESSL = max(abs(EulerErrSSL(2:nobs+1,:)))
RMSEESSL = sqrt(mean(EulerErrSSL(2:nobs+1,:).^2))
toc

%  current state linarization
tic;
[XCSL, ~, EulerErrCSL] = LinApp_CSL_Euler(@UnstFisc_dyn,param,X0',Z,...
    NN,logX,Eps,Phi);
kCSL = XCSL(:,1);
BCSL = XCSL(:,2);
yCSL = zeros(nobs+2,1);
iCSL = zeros(nobs+2,1);
cCSL = zeros(nobs+2,1);
rCSL = zeros(nobs+2,1);
wCSL = zeros(nobs+2,1);
tauCSL = zeros(nobs+2,1);
KCSL = zeros(nobs+2,1);
for t=1:nobs+1
    [yCSL(t+1),iCSL(t+1),cCSL(t+1),rCSL(t+1),wCSL(t+1),tauCSL(t+1)...
        ,KCSL(t+1)] = UnstFisc_defs(kCSL(t),BCSL(t),Z(t+1),kCSL(t+1),...
        BCSL(t+1),param);
end

AvgEECSL = mean(EulerErrCSL(2:nobs+1,:))
MaxAEECSL = max(abs(EulerErrCSL(2:nobs+1,:)))
RMSEECSL = sqrt(mean(EulerErrCSL(2:nobs+1,:).^2))
toc


% plot
plotdataB = [BSSL BCSL];
plotdatak = [kSSL(2:nobs+1) kCSL(2:nobs+1)];
plotdatatau = [tauSSL(2:nobs+1) tauCSL(2:nobs+1)];
plotdataBrat = [BSSL(2:nobs+1)./(kSSL(2:nobs+1) + BSSL(2:nobs+1)) ...
             BCSL(2:nobs+1)./(kCSL(2:nobs+1) + BCSL(2:nobs+1))];
UnstFiscal_Plots(plotdataB,plotdatak) 
PlotCompare2(plotdatatau) 
PlotCompare2(plotdataBrat) 