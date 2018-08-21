function [eulerrs, Kvec] = CSLStdInfHor(sim_type)

% This program simulates an otherwise standard RBC model with inelastically
% supplied labor.  No taxes

% The argument sim_type specifies whether we are simulating a 1 std dev
% shock to the model, a 2 std dev shock and starting at 1.5 times the ss
% or a random simulation.  The input must be either:
% '1stddev'
% '2stddev'
% 'randsim'
disp(sim_type);
global alpha beta delta gamma lbar sigz rhoz nx ny nz incr Zbar PP QQ UU
global Kzsp wzsp rzsp K0 czsp tauzsp T Kbar zzsp Yzsp discount

% set model parameters
alpha = .35;   %capital share in GDP .3
beta  = .96;  %discount factor .995
delta = .05; %depreciation rate 100% for exact solution to hold
lbar  = 1;     %labor per worker  1
gamma = 2;  %utility curvature
sigz = 0.0130058759370760;   %standard deviation for shocks .02
rhoz  = .95^4;    %autocorrelation of z series .9
zbar  = 0;
Zbar  = zbar;
rbar = (1 / beta) - 1 + delta;

% set starting state in terms of percent deviations from SS
Kstart = -.1;
zstart = 0;

% set program parameters
nmc = 1; %number of monte carlos
T = 1000; %sample size

%Solve for SS
guess = .03;
options = optimset('Display','iter','MaxFunEvals',100000000,'MaxIter', 200000 ...
				   ,'Algorithm','levenberg-marquardt','TolX',1.0e-8,'TolFun',1.0e-8);
Xbar = fsolve(@CSLStdInfHorss,guess,options);
guess = real(Xbar);
Xbar = fsolve(@CSLStdInfHorss,guess,options);
Kbar = Xbar;

temp = CSLStdInfHordefs(Xbar, zbar, Xbar);

Ybar = temp(1);
wbar = temp(2);
rbar = temp(3);
cbar = temp(4);

cguess = cbar;


% read in SS values of input vector
in = [Kbar; Kbar; Kbar; zbar; zbar];

% check SS solution by passing in to dynamic eqations and confirm =0
zz = CSLStdInfHordyn(in);

if abs(zz)>1e-8;
	disp('SS must be off because Euler Error is not within given tolerance')
end;

% get PP & QQ matrices for log-linearization about the steady state
nx=1; %number of endogenous state variables (x's)
ny=0; %number of jump variables (y's)
nz=1; %number of exogenous state variables (z's)
incr= 1e-6; % epsilon for computing numerical derivatives

% initial value of K
disp(sim_type)
if sim_type == '2stddev'; 
    K0 = 1.5*Kbar;

elseif sim_type == '1stddev' ;
    K0 = Kbar;

elseif sim_type == 'randsim' ;
    K0 = Kbar;
end

% Create a vector of parameters that will be passed into the calculate
% euler errors function
params = [beta, gamma, delta, 0, sigz, rhoz, alpha, Kbar]

% read in SS values of input vector
in = [Kbar; Kbar; Kbar; zbar; zbar];

eedyncsl = zeros(T,1) ;
% Toggle these comments for the simulation.
% load('simulationshocks')
% T = 150

% begin monte carlos
for m=1:nmc
    m ;
    %generate z series
    
    % Change up the eps depending on which simulation you're running.
    if sim_type == '2stddev'
    eps = zeros(T, 1);
    eps(2) = sigz * 2;
    
    elseif sim_type == '1stddev'
    eps = zeros(T, 1);
    eps(2) = sigz;
    
    elseif sim_type == 'randsim'
        load('simulationshocks/simulationshocks.mat');
        eps = simeps;
    end

    z = zeros(T,1);
    z(1) = rhoz * zstart ;
    for t=2:T
        z(t) = rhoz*z(t-1) + eps(t);
    end

    %generate appoximate time series using linearization about CS
    %calculating a new PP & QQ & UUeach period
    UUcase = 1;
    Kcs   = zeros(T,1);
    Ycs   = zeros(T,1);
    wcs   = zeros(T,1);
    rcs   = zeros(T,1);
    ccs   = zeros(T,1);
    PPcs  = zeros(T,1);
    QQcs  = zeros(T,1);
    UUcs  = zeros(T,1);
	
	
    Kcs(1) = K0;
	
    for t=1:T-1
        [Kcs(t+1),PPcs(t),QQcs(t),UUcs(t)] = CSL1get(Kcs(t),z(t));
        eedyncsl(t) = calceulerrs(Kcs(t), z(t), ...
                                  UUcs(t), params) ;
        temp     = CSLStdInfHordefs(Kcs(t),z(t),Kcs(t+1));
        Ycs(t)   = temp(1);
        wcs(t)   = temp(2);
        rcs(t)   = temp(3);
        ccs(t)   = temp(4);
	end
	
    [Kcslast,PPcs(T),QQcs(T),UUcs(T)] = CSL1get(Kcs(T),z(T));
    temp = CSLStdInfHordefs(Kcs(T),z(T),Kcslast);
    eedyncsl(t) = calceulerrs(Kcs(T), z(T), ...
                              UUcs(T), params) ;
    Ycs(T) = temp(1);
    wcs(T) = temp(2);
    rcs(T) = temp(3);
    ccs(T) = temp(4);
	
    %% find MAD of Kcs & Kss from Ktp
    %MADss = mean(abs(Kss-Ktp));
    %MADcs = mean(abs(Kcs-Ktp));
    %% add to MC average
    %MADssmc = (m-1)*MADssmc/m + MADss/m;
    %MADcsmc = (m-1)*MADcsmc/m + MADcs/m;
    %COcsmc = (m-1)*COcsmc/m + mean([PPcs QQcs UUcs])/m;

end
% TODO: RMSE.  Need to divide by n inside sqrt
RMSE = sqrt(sum(eedyncsl.^2)/T) ;
Maxee = max(abs(eedyncsl)) ;
Minee = min(abs(eedyncsl)) ;

eulerrs = [RMSE, Maxee, Minee] ;
Kvec = Kcs ;

% Toggle the comments on these depending on which type of simulation you
% are trying to run.

% save('kmatsim', 'Kcs')
% save('kmatirfcsl', 'Kcs')
% save('kmat2sdirfcsl', 'Kcs')

%Add these plots in later if we need them.
%[MADssmc MADcsmc]
%[COss; COcsmc; mean(COzsp)]
%
%% % plot time paths of capital stocks
%% figure;
%% AAdata = [Kss Kcs Ktp Kzsp Kbar*ones(T,1)];
%% plot(AAdata)
%% legend('Kss','Kcs','Ktp','Kzsp','Kbar','Location','EastOutside')
%
%% plot simulated time paths for K,Y,c,w,r & tau
%figure;
%subplot (2,3,1)
%data = [Kss Kcs];
%plot(data)
%legend('Kss','Kcs','Location','North')
%subplot (2,3,2)
%data = [Yss Ycs Ytp];
%plot(data)
%legend('Yss','Ycs','Ytp','Location','North')
%subplot (2,3,3)
%data = [css ccs ctp];
%plot(data)
%legend('css','ccs','ctp','Location','North')
%subplot (2,3,4)
%data = [wss wcs wtp];
%plot(data)
%legend('wss','wcs','wtp','Location','North')
%subplot (2,3,5)
%data = [rss rcs rtp] ;
%plot(data)
%legend('rss','rcs','rtp','Location','North')
%subplot (2,3,6)
%data = [tauss taucs tautp];
%plot(data)
%legend('tauss','taucs','tautp','Location','North')
%
%% plot time paths for linear coefficients
%figure;
%subplot (1,3,1)
%data = [PPss*ones(T,1) PPcs PPzsp];
%plot(data)
%subplot (1,3,2)
%data = [QQss*ones(T,1) QQcs QQzsp];
%plot(data)
%subplot (1,3,3)
%data = [UUss*ones(T,1) UUcs UUzsp];
%plot(data)