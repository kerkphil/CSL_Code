function [eulerrs, Kvec] = CSLClosedForm(sim_type)

% This program simulates an otherwise standard RBC model with inelastically
% supplied labor.  No taxes

% The argument sim_type specifies whether we are simulating a 1 std dev
% shock to the model, a 2 std dev shock and starting at 1.5 times the ss
% or a random simulation.  The input must be either:
% '1stddev'
% '2stddev'
% 'randsim'

%-------------------------------------------------------------------------%
% Compute the 1 std dev IRF of the Brock-Mirman model with known 
% analytical solution using Current State Linearization Approximation
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
starttime = clock ; % start timer

%-------------------------------------------------------------------------%
% Load Model Parameteres
%-------------------------------------------------------------------------%
load '../ParamSet/ClosedFormParams.mat' beta alpha L mu rho sigma T rss ...
     kss zvec zTrans ;
gamma       = 2
rhoz        = rho ;
sigz        = sigma ;
T           = 40 ;

% We use K(z)bar as the ss val in some of the following commands so we simply 
% rename the steady state value of k(z) that was imported as Kbar (zbar)
Kbar = kss
zbar = mu
zstart = 0

% read in SS values of input vector
in = [Kbar; Kbar; Kbar; zbar; zbar];

% check SS solution by passing in to dynamic eqations and confirm =0
zz = CSLClosedFormdyn(in);

if abs(zz)>1e-8;
	disp('SS must be off because Euler Error is not within given tolerance')
end;


%-------------------------------------------------------------------------%
% Set up the parameters necessary for the Current State Linearization code
%-------------------------------------------------------------------------%
nx=1; %number of endogenous state variables (x's)
ny=0; %number of jump variables (y's)
nz=1; %number of exogenous state variables (z's)
incr= 1e-6; % epsilon for computing numerical derivatives

% initial value of K
if sim_type == '2stddev'; 
    K0 = 1.5*Kbar;

elseif sim_type == '1stddev' ;
    K0 = Kbar;

elseif sim_type == 'randsim' ;
    K0 = Kbar;
end

% read in SS values of input vector
in = [Kbar; Kbar; Kbar; zbar; zbar];


%-------------------------------------------------------------------------%
% Here we begin using the Current State Linearization to calculate the
% values of only K in the following periods.  This is done using the eps
% and z vecs that were loaded previously.
%-------------------------------------------------------------------------%

% Create a vector of zeros to be filled with the capital values and set
% first value to K0 (in this case Kbar)
Kcs   = zeros(T,1);
Kcs(1) = K0;

% Change up the eps depending on which simulation you're running.
if sim_type == '2stddev'
eps = zeros(T, 1);
eps(2) = sigz * 2;

elseif sim_type == '1stddev'
eps = zeros(T, 1);
eps(2) = sigz;

elseif sim_type == 'randsim'
    load('simulationshocks/simulationshocks.mat');
    eps = simeps(1:40);
end

% Create a vector of parameters that will be passed into the calculate
% euler errors function
params = [beta, gamma, 1, 0, sigz, rhoz, alpha, Kbar]

z = zeros(T,1);
z(1) = rhoz * zstart ;
for t=2:T
    z(t) = rhoz*z(t-1) + eps(t);
end

eedyncsl = zeros(T,1) ;
for t=1:T
    [Kcs(t+1),PPcs(t),QQcs(t),UUcs(t)] = CSL1get(Kcs(t),z(t));
    eedyncsl(t) = calceulerrs(Kcs(t), z(t), ...
                                  UUcs(t), params) ;
end

%-------------------------------------------------------------------------%
% Calculate Euler errors from computed solution
%-------------------------------------------------------------------------%
% eulvec1  = 1 x T vector of Euler errors for periods t = 0,1,...T-1
% tind     = integer, index for time period
% zt       = scalar, current period z value
% kt       = scalar, current period capital value
% ktp1     = scalar, next period capital value
% wt       = scalar, current period real wage
% rt       = scalar, current period real return on capital
% ct       = scalar, current period consumption
% mut      = scalar, marginal utility of current period consumption
% zminind  = integer, index of value in discretized support of z to which
%            zt is closest
% zprobs   = 1 x znodes vector, probability weights for ztp1 tomorrow which
%            are a linear interpolation of the transition matrix row for
%            which zt is closest
% zuw      = scalar, weight on upper node of zvec for convex combination of
%            transition matrix rows
% zlw      = scalar, weight on lower node of zvec for convex combination of
%            transition matrix rows
% rtp1v    = 1 x znodes vector, vector of potential values for rtp1 next
%            period
% wtp1v    = 1 x znodes vector, vector of potential values for wtp1 next
%            period
% ktp2v    = 1 x znodes vector, vector of potential values for ktp2 next
%            period
% Emutp1   = scalar, discounted expected marginal utility next period
% dst1_rms = scalar, root mean squared euler error
% dst1_mab = scalar, maximum absolute euler error
%-------------------------------------------------------------------------%
% eulvec_cs = zeros(1,T) ;
% for tind = 1:T
%     zt = zirfvec(tind) ;
%     kt = Kcs(tind) ;
%     ktp1 = Kcs(tind+1) ;
%     wt = (1-alpha)*exp(zt)*(kt^alpha) ;
%     rt = alpha*exp(zt)*(kt^(alpha-1)) ;
%     ct = wt + rt*kt - ktp1 ;
%     mut = 1/ct ;
%     [~,zminind] = min((repmat(zt,[length(zvec),1]) - zvec).^2,[],1) ;
%     if zt <= zvec(1)
%         zprobs = zTrans(1,:) ;
%     elseif zt <= zvec(zminind)
%         zuw = (zt - zvec(zminind-1))/(zvec(zminind) - zvec(zminind-1)) ;
%         zlw = 1 - zuw ;
%         zprobs = zlw*zTrans(zminind-1,:) + zuw*zTrans(zminind,:) ;
%     elseif zt > zvec(zminind) && zt < zvec(length(zvec))
%         zlw = (zvec(zminind+1) - zt)/(zvec(zminind+1) - zvec(zminind)) ;
%         zuw = 1 - zlw ;
%         zprobs = zlw*zTrans(zminind,:) + zuw*zTrans(zminind+1,:) ;
%     elseif zt >= zvec(length(zvec))
%         zprobs = zTrans(length(zvec),:) ;
%     end
%     rtp1v = (alpha*ktp1^(alpha-1)).*exp(zvec') ;
%     wtp1v = ((1-alpha)*(ktp1^alpha)).*exp(zvec') ;
%     ktp2v = (alpha*beta*(ktp1^alpha)).*exp(zvec') ;
%     Emutp1 = sum(zprobs.*rtp1v./(wtp1v + ktp1*rtp1v - ktp2v)) ;
%     eulvec_cs(tind) = beta*Emutp1/mut - 1 ;
% end

% dstcs_rsm = sqrt(sum(eulvec_cs.^2)/T) ;
% dstcs_mab = max(abs(eulvec_cs)) ;
% dstcs_minab = min(abs(eulvec_cs)) ;

RMSE = sqrt(sum(eedyncsl.^2)/T) ;
Maxee = max(abs(eedyncsl)) ;
Minee = min(abs(eedyncsl)) ;

eulerrs = [RMSE, Maxee, Minee] ;
Kvec = Kcs

%-------------------------------------------------------------------------%
% Plot IRF
%-------------------------------------------------------------------------%
% figure(2)
% plot((0:1:T),Kcs)
% xlabel('Period (t)','FontSize', 12);
% ylabel('Aggregate Capital (K)','FontSize', 12);

%-------------------------------------------------------------------------%
% runtime = etime(clock,starttime) ; % end timer
% save ClosedForm_d1run.mat ;