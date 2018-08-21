%-------------------------------------------------------------------------%
% Compute the IRF of the Brock-Mirman model with known analytical solution
% using value function iteration
%-------------------------------------------------------------------------%
% This program calls the following .m and .mat file(s)
%   * ClosedFormParams.mat: loads in common parameters across solution
%        methods
%   * psisolve_vfi.m: this m-file solves for the policy function psi(k,z)
%        inside each iteration of the value function iteration
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/ClosedForm/VFI') ;

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
% zsize    = scalar, number of nodes in discrete support of z
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

%-------------------------------------------------------------------------%
% Set grid over state-space in k dimension
%-------------------------------------------------------------------------%
% ksize   = scalar, number of grid points in k dimension of state space
%           grid
% kmin    = scalar, minimum value of k
% kmax    = scalar, maximum value of k
% kvec    = ksize x 1 vector of grid points in k dimension
% kmat    = ksize x zsize matrix of kvec vector copied across zsize columns
% zmat    = ksize x zsize matrix of zvec vector copied down ksize rows
%-------------------------------------------------------------------------%
ksize   = 100 ;
kmin    = 0.90*kss ;
kmax    = 1.15*kss ;
kvec    = (kmin:(kmax-kmin)/(ksize-1):kmax)' ;
kmat    = repmat(kvec,[1,zsize]) ;
zmat    = repmat(zvec',[ksize,1]) ;

%-------------------------------------------------------------------------%
% Solve for equilibrium policy function and value function by VFI
%-------------------------------------------------------------------------%
% Vinit    = ksize x zsize matrix, initial guess for value function
% Vnew     = ksize x zsize matrix, new value function from Bellman operator
% Vdist    = scalar, norm distance measure between value functions
% Vtoler   = scalar, tolerance distance below which fixed point has
%            converged
% Viter    = integer index of iteration number of value function iteration
% Vmaxiter = integer, maximum number of iterations for VFI
% wmat     = ksize x zsize matrix, equilibrium real wage for state grid
% rmat     = ksize x zsize matrix, equilibrium real return for state grid
% optionsp = options for fmincon in value function iteration
% Vprime   = ksize x zsize matrix, value function next period
% psimat   = ksize x zsize matrix, policy function for kprime
% EVprmat  = ksize x zsize matrix, EV(k',z)
% kind     = integer index of node in support of k
% zind     = integer index of node in support of z
%-------------------------------------------------------------------------%
Vinit    = zeros(ksize,zsize) ;
Vnew     = Vinit ;
Vdist    = 10 ;
Vtoler   = 10^(-8) ;
Viter    = 0 ;
Vmaxiter = 500 ;
wmat     = (1-alpha)*exp(zmat).*(kmat.^alpha) ;
rmat     = alpha*exp(zmat).*(kmat.^(alpha-1)) ;
options = optimset('Display','off','MaxFunEvals',100000,'MaxIter',1000,...
                    'TolFun',1e-15,'Algorithm','active-set') ;
while Vdist > Vtoler && Viter < Vmaxiter
    Viter  = Viter + 1 ;
    Vprime = Vnew ;
    % Generate new policy function and value function
    psimat  = zeros(ksize,zsize) ;
    Vnew    = zeros(ksize,zsize) ;
    EVprmat = sum(repmat(reshape(Vprime,[ksize,1,zsize]),[1,zsize,1]).*...
              repmat(reshape(zTrans,[1,zsize,zsize]),[ksize,1,1]),3) ;
    for kind = 1:ksize
        for zind = 1:zsize
            %-------------------------------------------------------------%
            % kt      = scalar, current period value of k
            % wt      = scalar, current period value of w for particular k
            %           and z
            % rt      = scalar, current period value of r for particular k
            %           and z
            % EVprvec = ksize x 1 vector, the column of EVprmat EV(k',z)
            %           corresponding to zind
            % params  = 1 x 20 vector of parameters to pass in to fmincon
            %           solution for capital investment policy function
            % kprinit = scalar, initial guess for optimal kprime for
            %           particular (k,z)
            % kpr     = scalar, solution for optimal kprime for particular
            %           (k,z)
            % psimat  = ksize x zsize matrix, policy function for k'
            % ct      = scalar, current period value of c for particular k,
            %           z, and k'
            % uct     = scalar, current period utility of consumption ct
            %           for particular k, z, and k'
            % m       = scalar, slope coefficient of extrapolation line
            % b       = scalar, y-intercept of extrapolation line
            % EVpr    = scalar, interpolated value of E[V(k',z)]
            %-------------------------------------------------------------%
            kt      = kvec(kind) ;
            wt      = wmat(kind,zind) ;
            rt      = rmat(kind,zind) ;
            EVprvec = EVprmat(:,zind) ;
            params  = [beta kt wt rt] ;
            kprinit = (wt + rt*kt)/2 ;
            [kpr,~] = fmincon(@psisolve_vfi,kprinit,[],...
               [],[],[],10^(-10),wt+rt*kt-10^(-10),[],options,params,...
               EVprvec,kvec) ;
            psimat(kind,zind) = kpr ;
            ct      = wt + rt*kt - kpr ;
            uct     = log(ct) ;

            % Approximate grid for value function with linear interpolation
            if kpr <= kvec(1)
                % linear extrapolation from lowest two points
                m = (EVprvec(2) - EVprvec(1))/(kvec(2) - kvec(1)) ;
                b = -m*kvec(1) + EVprvec(1) ;
                EVpr = m*kpr + b ;
            elseif kpr > kvec(1) && kpr < kvec(length(kvec))
                % cubic spline interpolation between points on interior of
                % grid
                EVpr = interp1(kvec,EVprvec,kpr,'spline') ;
            elseif kpr >= kvec(length(kvec)) ;
                % linear extrapolation from highest two points
                m = (EVprvec(length(EVprvec)) - EVprvec(length(EVprvec)-1))/...
                    (kvec(length(kvec)) - kvec(length(kvec)-1)) ;
                b = -m*kvec(length(kvec)) + EVprvec(length(EVprvec)) ;
                EVpr = m*kpr + b ;
            end
            Vnew(kind,zind) = uct + beta*EVpr ;
        end
    end
    Vdist = sum((Vprime(:) - Vnew(:)).^2) ;
    disp(Viter) ;
    disp(Vdist) ;
end

%-------------------------------------------------------------------------%
% Generate two IRFs using VFI policy function k' = psi(k,z)
%-------------------------------------------------------------------------%
% kirfvec_vfi1 = 1 x (Tirf+1) vector, IRF1 path for capital K using state
%                at t=0 of (k_0,z_0) = (kss,mu)
% kirfvec_vfi2 = 1 x (Tirf+1) vector, IRF2 path for capital K using state
%                at t=0 of (k_0,z_0) = (1.2*kss,mu). Note: we have to force
%                k_1 = 1.1*kss
% eulvec_vfi1  = 1 x Tirf vector of Euler errors from IRF1
% eulvec_vfi2  = 1 x Tirf vector of Euler errors from IRF2
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
kirfvec_vfi1 = [kss,zeros(1,Tirf)] ;
kirfvec_vfi2 = [1.1*kss,zeros(1,Tirf)] ;
eulvec_vfi1  = zeros(1,Tirf) ;
eulvec_vfi2  = zeros(1,Tirf) ;
for tind = 1:Tirf
    zt_1 = zirfvec1(tind) ;
    zt_2 = zirfvec2(tind) ;
    kt_1 = kirfvec_vfi1(tind) ;
    kt_2 = kirfvec_vfi2(tind) ;

    if zt_1 <= zvec(1) && kt_1 <= kvec(1)
        ktp1_1 = psimat(1,1) ;
    elseif zt_1 >= zvec(zsize) && kt_1 >= kvec(ksize)
        ktp1_1 = psimat(ksize,zsize) ;
    elseif zt_1 <= zvec(1) && kt_1 >= kvec(ksize)
        ktp1_1 = psimat(ksize,1) ;
    elseif zt_1 >= zvec(zsize) && kt_1 <= kvec(1)
        ktp1_1 = psimat(1,zsize) ;
    elseif zt_1 <= zvec(1) && kt_1 > kvec(1) && kt_1 < kvec(ksize)
        ktp1_1 = interp1(kvec,psimat(:,1),kt_1,'spline') ;
    elseif zt_1 >= zvec(zsize) && kt_1 > kvec(1) && kt_1 < kvec(ksize)
        ktp1_1 = interp1(kvec,psimat(:,zsize),kt_1,'spline') ;
    elseif zt_1 > zvec(1) && zt_1 < zvec(zsize) && kt_1 <= kvec(1)
        ktp1_1 = interp1(zvec,psimat(1,:),zt_1,'spline') ;
    elseif zt_1 > zvec(1) && zt_1 < zvec(zsize) && kt_1 >= kvec(ksize)
        ktp1_1 = interp1(zvec,psimat(ksize,:),zt_1,'spline') ;
    elseif zt_1 > zvec(1) && zt_1 < zvec(zsize) && ...
      kt_1 > kvec(1) && kt_1 < kvec(ksize)
        ktp1_1 = interp2(zmat,kmat,psimat,zt_1,kt_1,'spline') ;
    end

    if zt_2 <= zvec(1) && kt_2 <= kvec(1)
        ktp1_2 = psimat(1,1) ;
    elseif zt_2 >= zvec(zsize) && kt_2 >= kvec(ksize)
        ktp1_2 = psimat(ksize,zsize) ;
    elseif zt_2 <= zvec(1) && kt_2 >= kvec(ksize)
        ktp1_2 = psimat(ksize,1) ;
    elseif zt_2 >= zvec(zsize) && kt_2 <= kvec(1)
        ktp1_2 = psimat(1,zsize) ;
    elseif zt_2 <= zvec(1) && kt_2 > kvec(1) && kt_2 < kvec(ksize)
        ktp1_2 = interp1(kvec,psimat(:,1),kt_2,'spline') ;
    elseif zt_2 >= zvec(zsize) && kt_2 > kvec(1) && kt_2 < kvec(ksize)
        ktp1_2 = interp1(kvec,psimat(:,zsize),kt_2,'spline') ;
    elseif zt_2 > zvec(1) && zt_2 < zvec(zsize) && kt_2 <= kvec(1)
        ktp1_2 = interp1(zvec,psimat(1,:),zt_2,'spline') ;
    elseif zt_2 > zvec(1) && zt_2 < zvec(zsize) && kt_2 >= kvec(ksize)
        ktp1_2 = interp1(zvec,psimat(ksize,:),zt_2,'spline') ;
    elseif zt_2 > zvec(1) && zt_2 < zvec(zsize) && ...
      kt_2 > kvec(1) && kt_2 < kvec(ksize)
        ktp1_2 = interp2(zmat,kmat,psimat,zt_2,kt_2,'spline') ;
    end

    kirfvec_vfi1(tind+1) = ktp1_1 ;
    kirfvec_vfi2(tind+1) = ktp1_2 ;
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
    elseif zt_1 >= zvec(length(zvec))
        zprobs_1 = zTrans(length(zvec),:) ;
    elseif zt_1 <= zvec(zminind_1) && zt_1 > zvec(1)
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
    end
    rtp1v_1 = (alpha*ktp1_1^(alpha-1)).*exp(zvec') ;
    wtp1v_1 = ((1-alpha)*(ktp1_1^alpha)).*exp(zvec') ;
    ktp2v_1 = (alpha*beta*(ktp1_1^alpha)).*exp(zvec') ;
    Emutp1_1 = sum(zprobs_1.*rtp1v_1./...
               (wtp1v_1 + ktp1_1*rtp1v_1 - ktp2v_1)) ;
    eulvec_vfi1(tind) = beta*Emutp1_1/mut_1 - 1 ;

    [~,zminind_2] = min((repmat(zt_2,[length(zvec),1]) - zvec).^2,[],1) ;
    if zt_2 <= zvec(1)
        zprobs_2 = zTrans(1,:) ;
    elseif zt_2 >= zvec(length(zvec))
        zprobs_2 = zTrans(length(zvec),:) ;
    elseif zt_2 <= zvec(zminind_2) && zt_2 > zvec(1)
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
    end
    rtp1v_2 = (alpha*ktp1_2^(alpha-1)).*exp(zvec') ;
    wtp1v_2 = ((1-alpha)*(ktp1_2^alpha)).*exp(zvec') ;
    ktp2v_2 = (alpha*beta*(ktp1_2^alpha)).*exp(zvec') ;
    Emutp1_2 = sum(zprobs_2.*rtp1v_2./...
               (wtp1v_2 + ktp1_2*rtp1v_2 - ktp2v_2)) ;
    eulvec_vfi2(tind) = beta*Emutp1_2/mut_2 - 1 ;
end

%-------------------------------------------------------------------------%
% Generate simulated time path for K using VFI policy function
% k' = psi(k,z)
%-------------------------------------------------------------------------%
% ksimvec_vfi3 = 1 x Tsim+1 vector, simulated path of capital K starting at
%                state (k_0,z_0) = (kss,mu)
% eulvec_vfi3  = 1 x Tsim vector of Euler errors from longer simulation
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
% Fit a polynomial k'=exp(a0+a1*k+a2*z+a3*k^2+a4*z^2+a5*k*z) to psimat in
% order to interpolate
ktp1vec = reshape(psimat,[ksize*zsize,1]) ;
ktvec = reshape(repmat(kvec,[1,zsize]),[ksize*zsize,1]) ;
ztvec = reshape(repmat(zvec',[ksize,1]),[ksize*zsize,1]) ;
xx = [ones(size(ktvec)),ktvec,ztvec,ktvec.^2,ztvec.^2,ktvec.*ztvec] ;
psicoefs = (xx'*xx)\(xx'*ktp1vec) ;

ksimvec_vfi3 = [kss,zeros(1,Tsim)] ;
eulvec_vfi3  = zeros(1,Tsim) ;
for tind = 1:Tsim
    zt_3 = zsimvec3(tind) ;
    kt_3 = ksimvec_vfi3(tind) ;
    if zt_3 <= zvec(1) && kt_3 <= kvec(1)
        ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
    elseif zt_3 >= zvec(zsize) && kt_3 >= kvec(ksize)
        ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
    elseif zt_3 <= zvec(1) && kt_3 >= kvec(ksize)
        ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
    elseif zt_3 >= zvec(zsize) && kt_3 <= kvec(1)
        ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
    elseif zt_3 <= zvec(1) && kt_3 > kvec(1) && kt_3 < kvec(ksize)
        ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
    elseif zt_3 >= zvec(zsize) && kt_3 > kvec(1) && kt_3 < kvec(ksize)
        ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
    elseif zt_3 > zvec(1) && zt_3 < zvec(zsize) && kt_3 <= kvec(1)
        ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
    elseif zt_3 > zvec(1) && zt_3 < zvec(zsize) && kt_3 >= kvec(ksize)
        ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
    elseif zt_3 > zvec(1) && zt_3 < zvec(zsize) && ...
      kt_3 > kvec(1) && kt_3 < kvec(ksize)
        ktp1_3 = interp2(zmat,kmat,psimat,zt_3,kt_3,'spline') ;
    end
    ksimvec_vfi3(tind+1) = ktp1_3 ;
    wt_3 = (1-alpha)*exp(zt_3)*(kt_3^alpha) ;
    rt_3 = alpha*exp(zt_3)*(kt_3^(alpha-1)) ;
    ct_3 = wt_3 + rt_3*kt_3 - ktp1_3 ;
    mut_3 = 1/ct_3 ;

    [~,zminind_3] = min((repmat(zt_3,[length(zvec),1]) - zvec).^2,[],1) ;
    if zt_3 <= zvec(1)
        zprobs_3 = zTrans(1,:) ;
    elseif zt_3 >= zvec(length(zvec))
        zprobs_3 = zTrans(length(zvec),:) ;
    elseif zt_3 <= zvec(zminind_3) && zt_3 > zvec(1)
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
    end
    rtp1v_3 = (alpha*ktp1_3^(alpha-1)).*exp(zvec') ;
    wtp1v_3 = ((1-alpha)*(ktp1_3^alpha)).*exp(zvec') ;
    ktp2v_3 = (alpha*beta*(ktp1_3^alpha)).*exp(zvec') ;
    Emutp1_3 = sum(zprobs_3.*rtp1v_3./...
               (wtp1v_3 + ktp1_3*rtp1v_3 - ktp2v_3)) ;
    eulvec_vfi3(tind) = beta*Emutp1_3/mut_3 - 1 ;
end

%-------------------------------------------------------------------------%
% Calculate two Euler error distance measures for the two IRFs and the one
% simulation
%-------------------------------------------------------------------------%
% eulerr_rms_vfi1 = scalar, root mean squared Euler error from IRF1
% eulerr_mab_vfi1 = scalar, maximum absolute Euler error from IRF1
% eulerr_rms_vfi2 = scalar, root mean squared Euler error from IRF2
% eulerr_mab_vfi2 = scalar, maximum absolute Euler error from IRF1
% eulerr_rms_vfi3 = scalar, root mean squared Euler error from simulation
% eulerr_mab_vfi3 = scalar, maximum absolute Euler error from simulation
% anlerr_rms_vfi1 = scalar, root mean squared deviation from analytical
%                   solution for IRF1
% anlerr_rms_vfi2 = scalar, root mean squared deviation from analytical
%                   solution for IRF2
% anlerr_rms_vfi3 = scalar, root mean squared deviation from analytical
%                   solution for simulation
%-------------------------------------------------------------------------%
eulerr_rms_vfi1 = sqrt(sum(eulvec_vfi1.^2)/Tirf) ;
eulerr_mab_vfi1 = max(abs(eulvec_vfi1)) ;
eulerr_rms_vfi2 = sqrt(sum(eulvec_vfi2.^2)/Tirf) ;
eulerr_mab_vfi2 = max(abs(eulvec_vfi2)) ;
eulerr_rms_vfi3 = sqrt(sum(eulvec_vfi3.^2)/Tsim) ;
eulerr_mab_vfi3 = max(abs(eulvec_vfi3)) ;
anlerr_rms_vfi1 = sqrt(sum((kirfvec_vfi1 - kirfvec_ans1).^2)/Tirf) ;
anlerr_rms_vfi2 = sqrt(sum((kirfvec_vfi2 - kirfvec_ans2).^2)/Tirf) ;
anlerr_rms_vfi3 = sqrt(sum((ksimvec_vfi3 - ksimvec_ans3).^2)/Tirf) ;

display([eulerr_rms_vfi1,eulerr_mab_vfi1,anlerr_rms_vfi1;...
    eulerr_rms_vfi2,eulerr_mab_vfi2,anlerr_rms_vfi2;...
    eulerr_rms_vfi3,eulerr_mab_vfi3,anlerr_rms_vfi3]) ;

%-------------------------------------------------------------------------%
% Plot IRFs and simulation
%-------------------------------------------------------------------------%
figure(1)
plot((0:1:Tirf),kirfvec_vfi1,(0:1:Tirf),kirfvec_vfi2) ;
xlabel('Period (t)','FontSize', 12) ;
ylabel('Aggregate Capital (K)','FontSize', 12) ;

figure(2)
plot((0:1:Tsim),ksimvec_vfi3) ;
xlabel('Period (t)','FontSize', 12) ;
ylabel('Aggregate Capital (K)','FontSize', 12) ;

%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save ClosedForm_vfi.mat ;
