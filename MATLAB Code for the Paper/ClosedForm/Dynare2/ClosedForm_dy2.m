%-------------------------------------------------------------------------%
% Compute dynare 2nd order solution and corresponding two IRFs and
% longer simulation for Brock-Mirman (1972) infinite horizon model for
% which the closed form solution is known
%-------------------------------------------------------------------------%
% This program calls the following .mod and .mat file(s)
%     * ClosedForm_d2.mod: runs the dynare .mod file for the 2nd order
%          approximation
%     * ClosedFormParams.mat: loads common parameter values into memory
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/ClosedForm/Dynare2') ;

starttime = clock ; % start timer

%-------------------------------------------------------------------------%
% Run Dynare 2nd order approximation and pull out policy function values
%-------------------------------------------------------------------------%
% policy function is in the form:
%   k(t+1) = a0 + a1*(k(t)-kss) + a2*(z(t-1)-mu) + a3*eps(t)
%          + a4*(k(t)-kss)^2 + a5*(k(t)-kss)*(z(t-1)-mu) + a6*(z(t-1)-mu)^2
%          + a7*eps(t)^2 + a8*(k(t)-kss)*eps(t) + a9*(z(t-1)-mu)*eps(t)
%      This can equivalently be written
%   k(t+1) = b0 + b1*k(t) + b2*z(t) + b3*k(t)^2 + b4*k(t)*z(t) + b5*z(t)^2
%      where b0 = a0 - a1*kss - a3*mu + a4*kss^2 + a8*kss*mu + a7*mu^2
%            b1 = a1 - 2*a4*kss - a8*mu, b2 = a3 - a8*kss - 2*a7*mu
%            b3 = a4, b4 = a8, b5 = a7
%            a2/a3 = a5/a8 = rho, a6/a7 = rho^2, a9/a7 = 2*rho
% a0 = scalar, k(t+1) policy function constant equals kss
% a1 = scalar, k(t+1) policy function coefficient on k(t) - kss
% a2 = scalar, k(t+1) policy function coefficient on z(t-1) - mu
% a3 = scalar, k(t+1) policy function coefficient on eps(t)
% a4 = scalar, k(t+1) policy function coefficient on (k(t) - kss)^2
% a5 = scalar, k(t+1) policy function coefficient on
%      (k(t) - kss)*(z(t-1) - mu)
% a6 = scalar, k(t+1) policy function coefficient on (z(t-1) - mu)^2
% a7 = scalar, k(t+1) policy function coefficient on eps(t)^2
% a8 = scalar, k(t+1) policy function coefficient on (k(t) - kss)*eps(t)
% a9 = scalar, k(t+1) policy function coefficient on (z(t-1) - mu)*eps(t)
%-------------------------------------------------------------------------%
dynare ClosedForm_d2.mod noclearall ;
a0 = oo_.dr.ys(2,1) + 0.5*oo_.dr.ghs2(2,1) ;
a1 = oo_.dr.ghx(2,1) ;
a2 = oo_.dr.ghx(2,2) ;
a3 = oo_.dr.ghu(2,1) ;
a4 = 0.5*oo_.dr.ghxx(2,1) ;
a5 = oo_.dr.ghxx(2,2) ; % is actually 0.5*C(1,2) + 0.5*C(1,3)
a6 = 0.5*oo_.dr.ghxx(2,4) ;
a7 = 0.5*oo_.dr.ghuu(2,1) ;
a8 = oo_.dr.ghxu(2,1) ;
a9 = oo_.dr.ghxu(2,2) ;

%-------------------------------------------------------------------------%
% Read in parameters and generate simplified policy function coefficients
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
% zsize    = scalar, number of nodes in discrete support of z
% b0       = scalar, simplified k(t+1) policy function constant
% b1       = scalar, simplified k(t+1) policy function coefficient on k(t)
% b2       = scalar, simplified k(t+1) policy function coefficient on z(t)
% b3       = scalar, simplified k(t+1) policy function coefficient on
%            k(t)^2
% b4       = scalar, simplified k(t+1) policy function coefficient on
%            k(t)z(t)
% b5       = scalar, simplified k(t+1) policy function coefficient on
%            z(t)^2
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
zsize = length(zvec) ;
b0 = a0 - a1*kss - a3*mu + a4*kss^2 + a8*kss*mu + a7*mu^2 ;
b1 = a1 - 2*a4*kss - a8*mu ;
b2 = a3 - a8*kss - 2*a7*mu ;
b3 = a4 ;
b4 = a8 ;
b5 = a7 ;

%-------------------------------------------------------------------------%
% Generate two IRFs using dynare 2nd order policy function:
%             k' = b0 + b1*k + b2*z + b3*k^2 + b4*k*z + b5*z^2
%-------------------------------------------------------------------------%
% kirfvec_dy21 = 1 x (Tirf+1) vector, IRF1 path for capital K using state
%                at t=0 of (k_0,z_0) = (kss,mu)
% kirfvec_dy22 = 1 x (Tirf+1) vector, IRF2 path for capital K using state
%                at t=0 of (k_0,z_0) = (1.2*kss,mu). Note: we have to force
%                k_1 = 1.1*kss
% eulvec_dy21  = 1 x Tirf vector of Euler errors from IRF1
% eulvec_dy22  = 1 x Tirf vector of Euler errors from IRF2
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
kirfvec_dy21 = [kss,zeros(1,Tirf)] ;
kirfvec_dy22 = [1.1*kss,zeros(1,Tirf)] ;
eulvec_dy21  = zeros(1,Tirf) ;
eulvec_dy22  = zeros(1,Tirf) ;
for tind = 1:Tirf
    zt_1 = zirfvec1(tind) ;
    zt_2 = zirfvec2(tind) ;
    kt_1 = kirfvec_dy21(tind) ;
    kt_2 = kirfvec_dy22(tind) ;
    ktp1_1 = b0 + b1*kt_1 + b2*zt_1 + b3*kt_1^2 + b4*kt_1*zt_1 + ...
             b5*zt_1^2 ;
    ktp1_2 = b0 + b1*kt_2 + b2*zt_2 + b3*kt_2^2 + b4*kt_2*zt_2 + ...
             b5*zt_2^2 ;
    kirfvec_dy21(tind+1) = ktp1_1 ;
    kirfvec_dy22(tind+1) = ktp1_2 ;
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
    ktp2v_1 = (b0 + b1*ktp1_1 * b3*ktp1_1^2)*ones(1,zsize) + b2*zvec' + ...
              b4*ktp1_1*zvec' + b5*(zvec'.^2) ;
    Emutp1_1 = sum(zprobs_1.*rtp1v_1./...
               (wtp1v_1 + ktp1_1*rtp1v_1 - ktp2v_1)) ;
    eulvec_dy21(tind) = beta*Emutp1_1/mut_1 - 1 ;

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
    ktp2v_2 = (b0 + b1*ktp1_2 * b3*ktp1_2^2)*ones(1,zsize) + b2*zvec' + ...
              b4*ktp1_2*zvec' + b5*(zvec'.^2) ;
    Emutp1_2 = sum(zprobs_2.*rtp1v_2./...
               (wtp1v_2 + ktp1_2*rtp1v_2 - ktp2v_2)) ;
    eulvec_dy22(tind) = beta*Emutp1_2/mut_2 - 1 ;
end

%-------------------------------------------------------------------------%
% Generate simulated time path for K using dynare 2nd order policy function
%             k' = b0 + b1*k + b2*z + b3*k^2 + b4*k*z + b5*z^2
%-------------------------------------------------------------------------%
% ksimvec_dy23 = 1 x Tsim+1 vector, simulated path of capital K starting at
%                state (k_0,z_0) = (kss,mu)
% eulvec_dy23  = 1 x Tsim vector of Euler errors from longer simulation
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
ksimvec_dy23 = [kss,zeros(1,Tsim)] ;
eulvec_dy23  = zeros(1,Tsim) ;
for tind = 1:Tsim
    zt_3 = zsimvec3(tind) ;
    kt_3 = ksimvec_dy23(tind) ;
    ktp1_3 = b0 + b1*kt_3 + b2*zt_3 + b3*kt_3^2 + b4*kt_3*zt_3 + ...
             b5*zt_3^2 ;
    ksimvec_dy23(tind+1) = ktp1_3 ;
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
    ktp2v_3 = (b0 + b1*ktp1_3 + b3*ktp1_3^2)*ones(1,zsize) + b2*zvec' + ...
              b4*ktp1_3*zvec' + b5*(zvec'.^2) ;
    Emutp1_3 = sum(zprobs_3.*rtp1v_3./...
               (wtp1v_3 + ktp1_3*rtp1v_3 - ktp2v_3)) ;
    eulvec_dy23(tind) = beta*Emutp1_3/mut_3 - 1 ;
end

%-------------------------------------------------------------------------%
% Calculate two Euler error distance measures for the two IRFs and the one
% simulation
%-------------------------------------------------------------------------%
% eulerr_rms_dy21 = scalar, root mean squared Euler error from IRF1
% eulerr_mab_dy21 = scalar, maximum absolute Euler error from IRF1
% eulerr_rms_dy22 = scalar, root mean squared Euler error from IRF2
% eulerr_mab_dy22 = scalar, maximum absolute Euler error from IRF1
% eulerr_rms_dy23 = scalar, root mean squared Euler error from simulation
% eulerr_mab_dy23 = scalar, maximum absolute Euler error from simulation
% anlerr_rms_dy21 = scalar, root mean squared deviation from analytical
%                   solution for IRF1
% anlerr_rms_dy22 = scalar, root mean squared deviation from analytical
%                   solution for IRF2
% anlerr_rms_dy23 = scalar, root mean squared deviation from analytical
%                   solution for simulation
%-------------------------------------------------------------------------%
eulerr_rms_dy21 = sqrt(sum(eulvec_dy21.^2)/Tirf) ;
eulerr_mab_dy21 = max(abs(eulvec_dy21)) ;
eulerr_rms_dy22 = sqrt(sum(eulvec_dy22.^2)/Tirf) ;
eulerr_mab_dy22 = max(abs(eulvec_dy22)) ;
eulerr_rms_dy23 = sqrt(sum(eulvec_dy23.^2)/Tsim) ;
eulerr_mab_dy23 = max(abs(eulvec_dy23)) ;
anlerr_rms_dy21 = sqrt(sum((kirfvec_dy21 - kirfvec_ans1).^2)/Tirf) ;
anlerr_rms_dy22 = sqrt(sum((kirfvec_dy22 - kirfvec_ans2).^2)/Tirf) ;
anlerr_rms_dy23 = sqrt(sum((ksimvec_dy23 - ksimvec_ans3).^2)/Tirf) ;

display([eulerr_rms_dy21,eulerr_mab_dy21,anlerr_rms_dy21;...
    eulerr_rms_dy22,eulerr_mab_dy22,anlerr_rms_dy22;...
    eulerr_rms_dy23,eulerr_mab_dy23,anlerr_rms_dy23]) ;

%-------------------------------------------------------------------------%
% Plot IRFs and simulation
%-------------------------------------------------------------------------%
close ;

figure(1)
plot((0:1:Tirf),kirfvec_dy21,(0:1:Tirf),kirfvec_dy22) ;
xlabel('Period (t)','FontSize', 12) ;
ylabel('Aggregate Capital (K)','FontSize', 12) ;

figure(2)
plot((0:1:Tsim),ksimvec_dy23) ;
xlabel('Period (t)','FontSize', 12) ;
ylabel('Aggregate Capital (K)','FontSize', 12) ;

%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save ClosedForm_dy2.mat ;
