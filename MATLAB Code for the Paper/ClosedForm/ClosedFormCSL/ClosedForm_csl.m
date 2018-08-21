%-------------------------------------------------------------------------%
% Compute closed form analytical solution and corresponding two IRFs and
% longer simulation for Brock-Mirman (1972) infinite horizon model for
% which the closed form solution is known
%-------------------------------------------------------------------------%
% This program calls the following .mat file(s)
%     * ClosedFormParams.mat: loads common parameter values into memory
%     * CSL1get.m : Calculates the value of k(t+1) given current state
%           values for k(t) and z(t) using a log-linear or simple linear
%           approximation of the transition function about the current
%           state k(t) and z(t).
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/ClosedForm/ClosedFormCSL') ;

starttime = clock ; % start timer

%-------------------------------------------------------------------------%
% Read in parameters
%-------------------------------------------------------------------------%
% beta     = per-period discount factor
% alpha    = capital share of income, production function parameter
% L        = individual endowment of inelastically supplied labor
% mu       = unconditional mean of AR(1) process for z
% rho      = per-period persistence in AR(1) process for z
% sigma    = per-period standard deviation of shock in z process
% Tirf     = scalar, number of periods to simulate in IRFs
% Tsim     = scalar, number of periods to simulate in simulation test
% kss      = scalar, certainty-equivalent steady-state capital stock
% zvec     = znodes x 1 vector, discretized support of z
% zTrans   = znodes x znodes matrix, Tauchen-Hussey Markov transition
%            matrix for z
% zirfvec1 = 1 x Tirf vector of productivity levels for IRF1
% zirfvec2 = 1 x Tirf vector of productivity levels for IRF2
% zsimvec3 = 1 x Tsim vector of productivity levels for the simulation
% kirfvec_ans1 = 1 x (Tirf+1) vector, IRF1 path for capital K using state
%                at t=0 of (k_0,z_0) = (kss,mu)
% kirfvec_ans2 = 1 x (Tirf+1) vector, IRF2 path for capital K using state
%                at t=0 of (k_0,z_0) = (1.2*kss,mu). Note: we have to force
%                k_1 = 1.1*kss
% ksimvec_ans3 = 1 x Tsim vector, simulate time series of the aggregate
%                capital stock using closed form solution
%-------------------------------------------------------------------------%
load '../ParamSet/ClosedFormParams.mat' beta alpha L mu rho sigma Tirf ...
     Tsim kss zvec zTrans zirfvec1 zirfvec2 zsimvec3 ;
load '../ClosedFormSol/ClosedForm_ans.mat' kirfvec_ans1 kirfvec_ans2 ...
     ksimvec_ans3 ;

%-------------------------------------------------------------------------%
% Set some parameters for Uhlig method
%-------------------------------------------------------------------------%
% nx      = integer, number of endogenous state variables (elements in X)
% ny      = integer, number of endogenous jump variables (elements in Y)
% nz      = integer, number of exogenous state variables (elements in Z)
% incr    = scalar, epsilon for computing numerical derivatives
% uparams = 1 x 9 vector of parameters to pass into function calculating U
%-------------------------------------------------------------------------%
nx      = 1 ;
ny      = 0 ;
nz      = 1 ;
incr    = 0.00001 ;
uparams = [alpha, beta, mu, rho, sigma, nx, ny, nz, incr] ;

%-------------------------------------------------------------------------%
% Solve for PP, QQ, UU, and k(t+1) using Uhlig method over a grid of
% possible combinations of k(t) and z(t), where the point about which the
% system is linearized is k(t) and z(t):
%       k(t+1) = k(t) + P * [k(t) - k(t)] + Q * [z(t) - z(t)] + U
%-------------------------------------------------------------------------%
% kuvecmin = scalar, lowest value k node
% kuvecmax = scalar, highest value k node
% kusize   = integer, number of nodes in k over which to solve model
% kuvec    = kusize x 1 vector, grid in k over which solve model
% kumat    = kusize x zusize matrix, kuvec copied over zusize columns, used
%            for interp2 function below
% zusize   = integer, number of nodes in z over which to solve model
% zuvec    = 1 x zusize vector, grid in z over which solve model
% zumat    = kusize x zusize matrix, zuvec copied down kusize rows, used
%            for interp2 function below
% ktp1grid = kusize x zusize matrix, CSL solution for ktp1 policy function
%            as function of the current state grid
% PPgrid   = kusize x zusize matrix, coefficient on [k(t)-k(t)] from CSL
%            solution as function of the current state grid
% QQgrid   = kusize x zusize matrix, coefficient on [z(t)-z(t)] from CSL
%            solution as function of the current state grid
% UUgrid   = kusize x zusize matrix, constant U from CSL solution as
%            function of the current state grid
% kind     = integer, index over k nodes
% zind     = integer, index over z nodes
%-------------------------------------------------------------------------%
% kuvecmin = min(ksimvec_ans3) ;
% kuvecmax = max(ksimvec_ans3) ;
% kusize   = 20 ;
% kuvec    = linspace(kuvecmin,kuvecmax,kusize)' ;
% zuvec    = zvec' ;
% zusize   = length(zuvec) ;
% kumat    = repmat(kuvec,[1,zusize]) ;
% zumat    = repmat(zuvec,[kusize,1]) ;
% ktp1grid = zeros(kusize,zusize) ;
% PPgrid   = zeros(kusize,zusize) ;
% QQgrid   = zeros(kusize,zusize) ;
% UUgrid   = zeros(kusize,zusize) ;
% for kind = 1:kusize
%     for zind = 1:zusize
%         kt = kuvec(kind) ;
%         zt = zuvec(zind) ;
%         [ktp1,PP,QQ,UU] = CSL1get(kt,zt,uparams) ;
%         ktp1grid(kind,zind) = ktp1 ;
%         PPgrid(kind,zind) = PP ;
%         QQgrid(kind,zind) = QQ ;
%         UUgrid(kind,zind) = UU ;
%     end
% end

%-------------------------------------------------------------------------%
% Generate two IRFs using linearization about the current state (CSL)
% solution method: k(t+1) = k(t) + P*[k(t)-k(t)] + Q*[z(t)-z(t)] + U
%-------------------------------------------------------------------------%
% kirfvec_csl1 = 1 x (Tirf+1) vector, IRF1 path for capital K using state
%                at t=0 of (k_0,z_0) = (kss,mu)
% kirfvec_csl2 = 1 x (Tirf+1) vector, IRF2 path for capital K using state
%                at t=0 of (k_0,z_0) = (1.2*kss,mu). Note: we have to force
%                k_1 = 1.1*kss
% eulvec_csl1  = 1 x Tirf vector of Euler errors from IRF1
% eulvec_csl2  = 1 x Tirf vector of Euler errors from IRF2
% tind         = scalar index of time period
% zt_1         = scalar, current period z value for IRF1
% zt_2         = scalar, current period z value for IRF2
% kt_1         = scalar, current period k value for IRF1
% kt_2         = scalar, current period k value for IRF2
% ktp1_1       = scalar, next period capital stock value for IRF1
% ktp1_2       = scalar, next period capital stock value for IRF2
% wt_1         = scalar, current period w value for IRF1
% wt_2         = scalar, current period w value for IRF2
% rt_1         = scalar, current period r value for IRF1
% rt_2         = scalar, current period r value for IRF2
% ct_1         = scalar, current period c value for IRF1
% ct_2         = scalar, current period c value for IRF2
% mut_1        = scalar, current period marginal utility for IRF1
% mut_2        = scalar, current period marginal utility for IRF2
% zminind_1    = integer, index of value in discretized support of z to
%                which zt_1 is closest
% zminind_2    = integer, index of value in discretized support of z to
%                which zt_2 is closest
% zprobs_1     = 1 x znodes vector, probability weights for ztp1_1 tomorrow
%                which are a linear interpolation of the transition matrix
%                row for which zt_1 is closest
% zprobs_2     = 1 x znodes vector, probability weights for ztp1_2 tomorrow
%                which are a linear interpolation of the transition matrix
%                row for which zt_2 is closest
% zuw_1        = scalar, weight on upper node of zvec for convex
%                combination of transition matrix rows for IRF1
% zuw_2        = scalar, weight on upper node of zvec for convex
%                combination of transition matrix rows for IRF2
% zlw_1        = scalar, weight on lower node of zvec for convex
%                combination of transition matrix rows for IRF1
% zlw_2        = scalar, weight on lower node of zvec for convex
%                combination of transition matrix rows for IRF2
% rtp1v_1      = 1 x znodes vector of potential values for rtp1 next period
%                for IRF1
% rtp1v_2      = 1 x znodes vector of potential values for rtp1 next period
%                for IRF2
% wtp1v_1      = 1 x znodes vector of potential values for wtp1 next period
%                for IRF1
% wtp1v_2      = 1 x znodes vector of potential values for wtp1 next period
%                for IRF2
% ktp2v_1      = 1 x znodes vector of potential values for ktp2 next period
%                for IRF1
% ktp2v_2      = 1 x znodes vector of potential values for ktp2 next period
%                for IRF2
% Emutp1_1     = scalar, discounted expected marginal utility next period
%                for IRF1
% Emutp1_2     = scalar, discounted expected marginal utility next period
%                for IRF2
%-------------------------------------------------------------------------%
kirfvec_csl1 = [kss,zeros(1,Tirf)] ;
kirfvec_csl2 = [1.1*kss,zeros(1,Tirf)] ;
%kirfvec_csl4 = [kss,zeros(1,Tirf)] ;
eulvec_csl1  = zeros(1,Tirf) ;
eulvec_csl2  = zeros(1,Tirf) ;
%eulvec_csl4  = zeros(1,Tirf) ;
for tind = 1:Tirf
    ktp2v_1 = zeros(size(zvec')) ;
    ktp2v_2 = zeros(size(zvec')) ;
    zt_1 = zirfvec1(tind) ;
    zt_2 = zirfvec2(tind) ;
    %zt_4 = zirfvec4(tind) ;
    kt_1 = kirfvec_csl1(tind) ;
    kt_2 = kirfvec_csl2(tind) ;
    %kt_4 = kirfvec_csl4(tind) ;
    [ktp1_1,PP,QQ,UU] = CSL1get(kt_1,zt_1,uparams) ;
    [ktp1_2,PP,QQ,UU] = CSL1get(kt_2,zt_2,uparams) ;
    %ktp1_4 = kt_4 + interp2(zumat,kumat,UUgrid,zt_4,kt_4) ;
    kirfvec_csl1(tind+1) = ktp1_1 ;
    kirfvec_csl2(tind+1) = ktp1_2 ;
    %kirfvec_csl4(tind+1) = ktp1_4 ;
    wt_1 = (1-alpha)*exp(zt_1)*(kt_1^alpha) ;
    wt_2 = (1-alpha)*exp(zt_2)*(kt_2^alpha) ;
    rt_1 = alpha*exp(zt_1)*(kt_1^(alpha-1)) ;
    rt_2 = alpha*exp(zt_2)*(kt_2^(alpha-1)) ;
    ct_1 = wt_1 + rt_1*kt_1 - ktp1_1 ;
    ct_2 = wt_2 + rt_2*kt_2 - ktp1_2 ;
    mut_1 = 1/ct_1 ;
    mut_2 = 1/ct_2 ;

    [~,zminind_1] = min((repmat(zt_1,[length(zvec),1]) - zvec).^2,[],1) ;
    if zt_1 <= zvec(1)
        zprobs_1 = zTrans(1,:) ;
    elseif zt_1 <= zvec(zminind_1)
        zuw_1    = (zt_1 - zvec(zminind_1-1))/...
                   (zvec(zminind_1) - zvec(zminind_1-1)) ;
        zlw_1    = 1 - zuw_1 ;
        zprobs_1 = zlw_1*zTrans(zminind_1-1,:) + ...
                   zuw_1*zTrans(zminind_1,:) ;
    elseif zt_1 > zvec(zminind_1) && zt_1 < zvec(length(zvec))
        zlw_1    = (zvec(zminind_1+1) - zt_1)/(zvec(zminind_1+1) - ...
                   zvec(zminind_1)) ;
        zuw_1    = 1 - zlw_1 ;
        zprobs_1 = zlw_1*zTrans(zminind_1,:) + ...
                   zuw_1*zTrans(zminind_1+1,:) ;
    elseif zt_1 >= zvec(length(zvec))
        zprobs_1 = zTrans(length(zvec),:) ;
    end
    rtp1v_1 = (alpha*ktp1_1^(alpha-1)).*exp(zvec') ;
    wtp1v_1 = ((1-alpha)*(ktp1_1^alpha)).*exp(zvec') ;
    for zind = 1:length(zvec)
        [ktp2v_1(zind),PP,QQ,UU] = CSL1get(ktp1_1,zvec(zind),uparams) ;
        [ktp2v_2(zind),PP,QQ,UU] = CSL1get(ktp1_2,zvec(zind),uparams) ;
    end
    Emutp1_1 = sum(zprobs_1.*rtp1v_1./...
               (wtp1v_1 + ktp1_1*rtp1v_1 - ktp2v_1)) ;
    eulvec_csl1(tind) = beta*Emutp1_1/mut_1 - 1 ;

    [~,zminind_2] = min((repmat(zt_2,[length(zvec),1]) - zvec).^2,[],1) ;
    if zt_2 <= zvec(1)
        zprobs_2 = zTrans(1,:) ;
    elseif zt_2 <= zvec(zminind_2)
        zuw_2    = (zt_2 - zvec(zminind_2-1))/...
                   (zvec(zminind_2) - zvec(zminind_2-1)) ;
        zlw_2    = 1 - zuw_2 ;
        zprobs_2 = zlw_2*zTrans(zminind_2-1,:) + ...
                   zuw_2*zTrans(zminind_2,:) ;
    elseif zt_2 > zvec(zminind_2) && zt_2 < zvec(length(zvec))
        zlw_2    = (zvec(zminind_2+1) - zt_2)/(zvec(zminind_2+1) - ...
                   zvec(zminind_2)) ;
        zuw_2    = 1 - zlw_2 ;
        zprobs_2 = zlw_2*zTrans(zminind_2,:) + ...
                   zuw_2*zTrans(zminind_2+1,:) ;
    elseif zt_2 >= zvec(length(zvec))
        zprobs_2 = zTrans(length(zvec),:) ;
    end
    rtp1v_2 = (alpha*ktp1_2^(alpha-1)).*exp(zvec') ;
    wtp1v_2 = ((1-alpha)*(ktp1_2^alpha)).*exp(zvec') ;
    Emutp1_2 = sum(zprobs_2.*rtp1v_2./...
               (wtp1v_2 + ktp1_2*rtp1v_2 - ktp2v_2)) ;
    eulvec_csl2(tind) = beta*Emutp1_2/mut_2 - 1 ;
end

%-------------------------------------------------------------------------%
% Generate simulated time path for K using linearization about the current
% state (CSL) solution method:
%            k(t+1) = k(t) + P*[k(t)-k(t)] + Q*[z(t)-z(t)] + U
%-------------------------------------------------------------------------%
% ksimvec_csl3 = 1 x Tsim+1 vector, simulated path of capital K starting at
%                state (k_0,z_0) = (kss,mu)
% eulvec_csl3  = 1 x Tsim vector of Euler errors from longer simulation
% zt_3         = scalar, current period z value for simulation
% kt_3         = scalar, current period k value for simulation
% ktp1_3       = scalar, next period capital stock value for simulation
% wt_3         = scalar, current period w value for simulation
% rt_3         = scalar, current period r value for simulation
% ct_3         = scalar, current period c value for simulation
% mut_3        = scalar, current period marginal utility for simulation
% zminind_3    = integer, index of value in discretized support of z to
%                which zt_3 is closest
% zprobs_3     = 1 x znodes vector, probability weights for ztp1_3 tomorrow
%                which are a linear interpolation of the transition matrix
%                row for which zt_3 is closest
% zuw_3        = scalar, weight on upper node of zvec for convex
%                combination of transition matrix rows for simulation
% zlw_3        = scalar, weight on lower node of zvec for convex
%                combination of transition matrix rows for simulation
% rtp1v_3      = 1 x znodes vector of potential values for rtp1 next period
%                for simulation
% wtp1v_3      = 1 x znodes vector of potential values for wtp1 next period
%                for simulation
% ktp2v_3      = 1 x znodes vector of potential values for ktp2 next period
%                for simulation
% Emutp1_3     = scalar, discounted expected marginal utility next period
%                for simulation
%-------------------------------------------------------------------------%
ksimvec_csl3 = [kss,zeros(1,Tsim)] ;
eulvec_csl3  = zeros(1,Tsim) ;
for tind = 1:Tsim
    ktp2v_3 = zeros(size(zvec')) ;
    zt_3 = zsimvec3(tind) ;
    kt_3 = ksimvec_csl3(tind) ;
    [ktp1_3,PP,QQ,UU] = CSL1get(kt_3,zt_3,uparams) ;
%     if zt_3 <= zuvec(1) && kt_3 <= kuvec(1)
%         ktp1_3 = ktp1grid(1,1) ;
%     elseif zt_3 >= zuvec(zusize) && kt_3 >= kuvec(kusize)
%         ktp1_3 = ktp1grid(kusize,zusize) ;
%     elseif zt_3 <= zuvec(1) && kt_3 >= kuvec(kusize)
%         ktp1_3 = ktp1grid(kusize,1) ;
%     elseif zt_3 >= zuvec(zusize) && kt_3 <= kuvec(1)
%         ktp1_3 = ktp1grid(1,zsize) ;
%     elseif zt_3 <= zuvec(1) && kt_3 > kuvec(1) && kt_3 < kuvec(kusize)
%         ktp1_3 = interp1(kuvec,ktp1grid(:,1),kt_3,'spline') ;
%     elseif zt_3 >= zuvec(zusize) && kt_3 > kuvec(1) && kt_3 < kuvec(kusize)
%         ktp1_3 = interp1(kuvec,ktp1grid(:,zusize),kt_3,'spline') ;
%     elseif zt_3 > zuvec(1) && zt_3 < zuvec(zusize) && kt_3 <= kuvec(1)
%         ktp1_3 = interp1(zuvec,ktp1grid(1,:),zt_3,'spline') ;
%     elseif zt_3 > zuvec(1) && zt_3 < zuvec(zusize) && kt_3 >= kuvec(kusize)
%         ktp1_3 = interp1(zuvec,ktp1grid(kusize,:),zt_3,'spline') ;
%     elseif zt_3 > zuvec(1) && zt_3 < zuvec(zusize) && ...
%       kt_3 > kuvec(1) && kt_3 < kuvec(kusize)
%         ktp1_3 = kt_3 + interp2(zumat,kumat,UUgrid,zt_3,kt_3) ;
%     end
    ksimvec_csl3(tind+1) = ktp1_3 ;
    wt_3 = (1-alpha)*exp(zt_3)*(kt_3^alpha) ;
    rt_3 = alpha*exp(zt_3)*(kt_3^(alpha-1)) ;
    ct_3 = wt_3 + rt_3*kt_3 - ktp1_3 ;
    mut_3 = 1/ct_3 ;

    [~,zminind_3] = min((repmat(zt_3,[length(zvec),1]) - zvec).^2,[],1) ;
    if zt_3 <= zvec(1)
        zprobs_3 = zTrans(1,:) ;
    elseif zt_3 <= zvec(zminind_3)
        zuw_3    = (zt_3 - zvec(zminind_3-1))/...
                   (zvec(zminind_3) - zvec(zminind_3-1)) ;
        zlw_3    = 1 - zuw_3 ;
        zprobs_3 = zlw_3*zTrans(zminind_3-1,:) + ...
                   zuw_3*zTrans(zminind_3,:) ;
    elseif zt_3 > zvec(zminind_3) && zt_3 < zvec(length(zvec))
        zlw_3    = (zvec(zminind_3+1) - zt_3)/(zvec(zminind_3+1) - ...
                   zvec(zminind_3)) ;
        zuw_3    = 1 - zlw_3 ;
        zprobs_3 = zlw_3*zTrans(zminind_3,:) + ...
                   zuw_3*zTrans(zminind_3+1,:) ;
    elseif zt_3 >= zvec(length(zvec))
        zprobs_3 = zTrans(length(zvec),:) ;
    end
    rtp1v_3 = (alpha*ktp1_3^(alpha-1)).*exp(zvec') ;
    wtp1v_3 = ((1-alpha)*(ktp1_3^alpha)).*exp(zvec') ;
    for zind = 1:length(zvec)
        [ktp2v_3(zind),PP,QQ,UU] = CSL1get(ktp1_3,zvec(zind),uparams) ;
    end
%     ktp1v_3 = ktp1_3*ones(size(zvec')) ;
%     ktp2v_3 = ktp1v_3 + interp2(zumat,kumat,UUgrid,zvec',ktp1v_3) ;
    Emutp1_3 = sum(zprobs_3.*rtp1v_3./...
               (wtp1v_3 + ktp1_3*rtp1v_3 - ktp2v_3)) ;
    eulvec_csl3(tind) = beta*Emutp1_3/mut_3 - 1 ;
end


%-------------------------------------------------------------------------%
% Calculate two Euler error distance measures for the two IRFs and the one
% simulation
%-------------------------------------------------------------------------%
% eulerr_rms_csl1 = scalar, root mean squared Euler error from IRF1
% eulerr_mab_csl1 = scalar, maximum absolute Euler error from IRF1
% eulerr_rms_csl2 = scalar, root mean squared Euler error from IRF2
% eulerr_mab_csl2 = scalar, maximum absolute Euler error from IRF1
% eulerr_rms_csl3 = scalar, root mean squared Euler error from simulation
% eulerr_mab_csl3 = scalar, maximum absolute Euler error from simulation
% anlerr_rms_csl1 = scalar, root mean squared deviation from analytical
%                   solution for IRF1
% anlerr_rms_csl2 = scalar, root mean squared deviation from analytical
%                   solution for IRF2
% anlerr_rms_csl3 = scalar, root mean squared deviation from analytical
%                   solution for simulation
%-------------------------------------------------------------------------%
eulerr_rms_csl1 = sqrt(sum(eulvec_csl1.^2)/Tirf) ;
eulerr_mab_csl1 = max(abs(eulvec_csl1)) ;
eulerr_rms_csl2 = sqrt(sum(eulvec_csl2.^2)/Tirf) ;
eulerr_mab_csl2 = max(abs(eulvec_csl2)) ;
eulerr_rms_csl3 = sqrt(sum(eulvec_csl3.^2)/Tsim) ;
eulerr_mab_csl3 = max(abs(eulvec_csl3)) ;
anlerr_rms_csl1 = sqrt(sum((kirfvec_csl1 - kirfvec_ans1).^2)/Tirf) ;
anlerr_rms_csl2 = sqrt(sum((kirfvec_csl2 - kirfvec_ans2).^2)/Tirf) ;
anlerr_rms_csl3 = sqrt(sum((ksimvec_csl3 - ksimvec_ans3).^2)/Tirf) ;

display([eulerr_rms_csl1,eulerr_mab_csl1 anlerr_rms_csl1;...
    eulerr_rms_csl2,eulerr_mab_csl2 anlerr_rms_csl2;...
    eulerr_rms_csl3,eulerr_mab_csl3 anlerr_rms_csl3]) ;

%-------------------------------------------------------------------------%
% Plot IRFs and simulation
%-------------------------------------------------------------------------%
figure(1)
plot((0:1:Tirf),kirfvec_csl1,(0:1:Tirf),kirfvec_csl2) ;
xlabel('Period (t)','FontSize', 12) ;
ylabel('Aggregate Capital (K)','FontSize', 12) ;

figure(2)
plot((0:1:Tsim),ksimvec_csl3) ;
xlabel('Period (t)','FontSize', 12) ;
ylabel('Aggregate Capital (K)','FontSize', 12) ;

%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save ClosedForm_csl.mat ;
