

% set model parameters
alpha = .35;   %capital share in GDP .3
beta  = .96;  %discount factor .995
delta = .05; %depreciation rate 100% for exact solution to hold
lbar  = 1;     %labor per worker  1
gamma = 2;  %utility curvature
sigz = 0.0130058759370760;   %standard deviation for shocks .02
rhoz  = .95^4;    %autocorrelation of z series .9

% Create a vector of parameters that will be passed
param = [alpha, beta, delta, gamma]

% set numerical parameters
xi = .95;  %convex weight on PP in update
gscale = .0000001; %scalar to multiply eye(nx) and get initital guess for PP 
Kscale = 1; %scaler to mutiply Kbar and get K0
Zscale = 0; %scaler to add to Zbar and get Z0
meps = 2.2E-16;  % Machine epsilon for double precision 2.2E-16
nx=1; %number of endogenous state variables (x's)
ny=0; %number of jump variables (y's)
nz=1; %number of exogenous state variables (z's)

% other values needed for calculation
zbar  = 0;
Zbar  = zbar;
NN = rhoz;

%Solve for SS
af = @(x) Mod1ss(x,param);
guess = .03;
options = optimset('Display','iter','MaxFunEvals',100000000,...
          'MaxIter', 200000,'Algorithm','levenberg-marquardt',...
          'TolX',1.0e-8,'TolFun',1.0e-8);
Xbar = fsolve(af,guess,options);
guess = real(Xbar);
Xbar = fsolve(af,guess,options);
Kbar = Xbar;

% read in SS values of input vector
in = [Kbar; Kbar; Kbar; zbar; zbar];

% check SS solution by passing in to dynamic eqations and confirm =0
zz = Mod1dyn(in,param);

if abs(zz)>1e-8;
    zz
	disp('SS must be off because Euler Error is not within given tolerance')
end;

% get PP & QQ matrices for log-linearization about the steady state

% linearizarion value of K
K0 = Kbar*Kscale;
Z0 = Zbar+Zscale;
% read in values of input vector
in = [K0; K0; K0; Z0; Z0];
% eps = max(sqrt(meps)*max(in,1000*sqrt(meps)*ones(5,1))) %value for derivatives
eps = 1e-7

% Find coefficients about K0
[AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WW, TT] = ...
    LinApp_Deriv(@Mod1dyn,param,in,nx,ny,nz,0);
[PP, QQ, UU, RR, SS, VV] = ...
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WW,TT,NN,0,0)

% read in plus values of input vector
inp = [K0+eps; K0+eps; K0+eps; Z0; Z0];

% Find coefficients about K0 plus
[AAp, BBp, CCp, DDp, FFp, GGp, HHp, JJp, KKp, LLp, MMp, WWp, TTp] = ...
    LinApp_Deriv(@Mod1dyn,param,inp,nx,ny,nz,0);
[PPp, QQp, UUp, RRp, SSp, VVp] = ...
    LinApp_Solve(AAp,BBp,CCp,DDp,FFp,GGp,HHp,JJp,KKp,LLp,MMp,WWp,TTp,...
    NN,0,0);

% read in minus SS values of input vector
inm = [K0-eps; K0-eps; K0-eps; Z0; Z0];

% Find coefficients about K0 minus
[AAm, BBm, CCm, DDm, FFm, GGm, HHm, JJm, KKm, LLm, MMm, WWm, TTm] = ...
    LinApp_Deriv(@Mod1dyn,param,inm,nx,ny,nz,0);
[PPm, QQm, UUm, RRm, SSm, VVm] = ...
    LinApp_Solve(AAm,BBm,CCm,DDm,FFm,GGm,HHm,JJm,KKm,LLm,MMm,WWm,TTm,...
    NN,0,0);
% calculate new value for PP
PPnew = (UUp-UUm)/(2*eps) + 1

% Starting guess for PP
PPnew = gscale*eye(nx);
% Initialize convergence statistic and count
conv = 1;
count = 0;

while conv>1e-8
    count = count + 1;
    PP = xi*PP + (1-xi)*PPnew;
    QQp = -(MMp+LLp*NN)/(FFp*NN+FFp*PP+GGp);
    UUp = -(FFp*PP+FFp+GGp) \ (TTp+(FFp*QQp+LLp)*(NN*Z0-Z0));
    QQm = -(MMm+LLm*NN)/(FFm*NN+FFm*PP+GGm);
    UUm = -(FFm*PP+FFm+GGm) \ (TTm+(FFm*QQm+LLm)*(NN*Z0-Z0));
    PPnew = (UUp-UUm)/(2*eps) + 1;
    conv = max(max(abs(PPnew - PP)));
end
count
PP = xi*PP + (1-xi)*PPnew
QQ = -(MM+LL*NN)/(FF*NN+FF*PP+GG)
UU = -(FF*PP+FF+GG) \ (TT+(FF*QQ+LL)*(NN*Z0-Z0))