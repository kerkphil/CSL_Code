%-------------------------------------------------------------------------%
% Solve the growth model with multiple steady-state marginal tax rates
% using policy function iteration optimizing with the Euler equation
%-------------------------------------------------------------------------%
% This program calls the following .m and .mat file(s)
%   * unstab_par.mat: loads in common parameters across solution methods
%   * unstabEuler_vfi.m: this m-file solves for the policy function
%        psi(k,z) using a zero finder on the Euler equation inside each
%        iteration of the value function iteration
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/unstab/VFI') ;

starttime = clock ; % start timer

%-------------------------------------------------------------------------%
% Read in parameters
%-------------------------------------------------------------------------%
% beta     = per-period discount factor
% alpha    = capital share of income, production function parameter
% delta    = depreciation rate of capital
% gamma    = coefficient of relative risk aversion
% L        = individual endowment of inelastically supplied labor
% mu       = unconditional mean of AR(1) process for z
% rho      = per-period persistence in AR(1) process for z
% sigma    = per-period standard deviation of shock in z process
% Bmin     = scalar, minimum cutoff value of government debt, below which
%            the marginal tax rate is 0.99
% Blow     = scalar, low cutoff value of government debt. The marginal tax
%            rate between Bmin and Blow is a linear function between tauss
%            and 0.99
% Bhigh    = scalar, high cutoff value of government savings. The marginal
%            tax rate between Bhigh and Bmax is a linear function between
%            tauss and 0
% Bmax     = scalar, maximum cutoff value of government savings, above
%            which the marginal tax rate is zero
% d        = scalar, fixed lump-sum transfer
% tauss    = scalar, certainty-equivalent unstable steady-state marginal
%            tax rate
% Tirf     = scalar, number of periods to simulate in IRFs
% Tsim     = scalar, number of periods to simulate in simulation test
% kss      = scalar, certainty-equivalent unstable steady-state captl stock
% zvec     = znodes x 1 vector, discretized support of z
% zTrans   = znodes x znodes matrix, Tauchen-Hussey Markov transition
%            matrix for z
% zirfvec1 = 1 x Tirf vector of productivity levels for IRF1
% zirfvec2 = 1 x Tirf vector of productivity levels for IRF2
% zsimvec3 = 1 x Tsim vector of productivity levels for the simulation
% zsize    = scalar, number of nodes in discrete support of z
%-------------------------------------------------------------------------%
load '../ParamSet/unstab_par.mat' beta alpha delta gamma L mu rho ...
     sigma Bmin Blow Bhigh Bmax d tauss Tirf Tsim kss zvec zTrans ...
     zirfvec1 zirfvec2 zsimvec3 ;
zsize = length(zvec) ;

%-------------------------------------------------------------------------%
% Set grid over state-space in k dimension and B dimension
%-------------------------------------------------------------------------%
% ksize   = scalar, number of grid points in k dimension of state grid
% Bsize   = scalar, number of grid points in B dimension of state grid
% kmin    = scalar, minimum value of k
% kmax    = scalar, maximum value of k
% kvec    = ksize x 1 vector of grid points in k dimension
% Bvec1   = 1 x tausize1 vector, government budget sizes for budget levels
%           between and including Bmin and Blow
% Bvec2   = 1 x tausize2 vector, government budget sizes for budget levels
%           between and not including Blow and Bhigh
% Bvec3   = 1 x tausize3 vector, government budget sizes for budget levels
%           between and including Bhigh and Bmax
% Bvec    = Bsize x 1 vector of grid points in B dimension
% tauvec1 = 1 x tausize1 vector, marginal tax rates for budget levels
%           between and including Bmin and Blow
% tauvec2 = 1 x tausize2 vector, marginal tax rates (constant tauss) for
%           budget levels between and not including Blow and Bhigh
% tauvec3 = 1 x tausize3 vector, marginal tax rates for budget levels
%           between and including Bhigh and Bmax
% tauvec  = Bsize x 1 vector, marginal tax rates corresponding to Bvec
% karr    = ksize x Bsize x zsize array of kvec vector copied across Bsize
%           columns and zsize deep
% Barr    = ksize x Bsize x zsize array of Bvec vector copied across ksize
%           rows and zsize deep
% zarr    = ksize x Bsize x zsize array of zvec vector copied down ksize
%           rows and Bsize columns
% tauarr  = ksize x Bsize x zsize array of tauvec vector copied across
%           ksize rows and zsize deep
%-------------------------------------------------------------------------%
ksize = 200 ;
Bsize = 200 ;
kmin = 0.8*kss ;
kmax = 1.2*kss ;
kvec = (kmin:(kmax-kmin)/(ksize-1):kmax)' ;
Bvec1 = linspace(Bmin,Blow,ceil(0.35*Bsize)) ;
Bvec2 = linspace(Blow+0.5*(Blow-Bmin)/ceil(0.35*Bsize),...
        Bhigh-0.5*(Blow-Bmin)/ceil(0.35*Bsize),Bsize-2*ceil(0.35*Bsize)) ;
Bvec3 = linspace(Bhigh,Bmax,ceil(0.35*Bsize)) ;
Bvec = [Bvec1,Bvec2,Bvec3]' ;
tauvec1 = ((0.99-tauss)/(Bmin-Blow))*Bvec1 + ...
          (tauss - Blow*((0.99-tauss)/(Bmin-Blow)))*ones(size(Bvec1)) ;
tauvec2 = tauss*ones(size(Bvec2)) ;
tauvec3 = (tauss/(Bhigh-Bmax))*Bvec3 + ...
          (tauss - Bhigh*(tauss/(Bhigh-Bmax)))*ones(size(Bvec3)) ;
tauvec = [tauvec1,tauvec2,tauvec3]' ;
karr = repmat(kvec,[1,Bsize,zsize]) ;
Barr = repmat(Bvec',[ksize,1,zsize]) ;
zarr = repmat(reshape(zvec,[1,1,zsize]),[ksize,Bsize,1]) ;
tauarr = repmat(tauvec',[ksize,1,zsize]) ;

%-------------------------------------------------------------------------%
% Solve for equilibrium policy function and value function by VFI
%-------------------------------------------------------------------------%
% psidist    = scalar, norm distance measure between value functions
% psitoler   = scalar, tolerance distance below which fixed point has
%              converged
% psiiter    = integer index of iteration number of VFI
% psimaxiter = integer, maximum number of iterations for VFI
% Karr       = ksize x Bsize x zsize array, equilibrium aggregate capital
%              stock for state grid
% warr       = ksize x Bsize x zsize array, equilibrium real wage for state
%              grid
% rarr       = ksize x Bsize x zsize array, equilibrium real return for
%              state grid
% onearr     = ksize x Bsize x zsize array of ones
% psi_init   = ksize x Bsize x zsize array, last period optimal policy
%              function (zeros): psi_T(k_T,B_T,z_T)
% psiarr     = ksize x Bsize x zsize array, policy function for kprime
% options    = options for fmincon in value function iteration
% psiprime   = ksize x Bsize x zsize array, policy function for k'':
%              k'' = psi'(k',B',z')
% cprime     = ksize x Bsize x zsize array, optimal consumpt. next period:
%              c'(k',B',z',psiprime)
% muprime    = ksize xBsize x zsize array, optimal marginal utility of
%              consumption next period:
%              [1 + r(k',B',z') - delta]*u'(c'(k',B',z',psiprime))
% Emuprarr   = ksize x Bsize x zsize array, expected value of muprime:
%              E[mu(k',B',z)]
% kind     = integer index of node in support of k
% zind     = integer index of node in support of z
%-------------------------------------------------------------------------%
psidist    = 10 ;
psitoler   = 10^(-10) ;
psiiter    = 0 ;
psimaxiter = 500 ;
Karr = karr + Barr ;
warr = (1-alpha)*exp(zarr).*(Karr.^alpha) ;
rarr = alpha*exp(zarr).*(Karr.^(alpha-1)) ;
onearr = ones(size(rarr)) ;
psi_init = zeros(ksize,Bsize,zsize) ;
psiarr = psi_init ;

options = optimset('Display','off','MaxFunEvals',100000,'MaxIter',1000,...
                   'TolFun',1e-15) ;

while psidist > psitoler && psiiter < psimaxiter
    psiiter  = psiiter + 1 ;
    psiprime = psiarr ;
    cprime = (onearr - tauarr).*warr + ((1 - delta)*onearr + rarr).*karr...
             - psiprime + d*onearr ;
    muprime = (onearr + rarr - delta*onearr).*(cprime.^(-gamma)) ;

    Emuprarr = sum( ...
        repmat(reshape(muprime,[ksize,Bsize,1,zsize]),[1,1,zsize,1]).* ...
        repmat(reshape(zTrans,[1,1,zsize,zsize]),[ksize,Bsize,1,1]),4) ;
    psiarr = zeros(ksize,Bsize,zsize) ;
           
    for kind = 1:ksize
        for Bind = 1:Bsize
            for zind = 1:zsize
                %-------------------------------------------------------------%
                % kt       = scalar, current period value of k
                % Bt       = scalar, current period value of B
                % taut     = scalar, current period value of tau
                % wt       = scalar, current period value of w for
                %            particular k, B, and z
                % rt       = scalar, current period value of r for
                %            particular k, B, and z
                % Bpr      = scalar, next period value of B for particular
                %            k, B, and z
                % Emuprvec = ksize x 1 vector, the column of Emuprmat
                %            E[mu(k',B',z)] corresponding to B' and zind
                % params  = 1 x 8 vector of parameters to pass in to fsolve
                % kprinit = scalar, initial guess for optimal kprime for
                %           particular (k,z)
                % kprmin  = scalar, minimum possible value for k'
                % kprmax  = scalar, maximum possible value for k'
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
                kt = kvec(kind) ;
                Bt = Bvec(Bind) ;
                taut = tauvec(Bind) ;
                wt = warr(kind,Bind,zind) ;
                rt = rarr(kind,Bind,zind) ;
                Bpr = taut*wt - d + (1+rt-delta)*Bt ;
                if kt + Bt <= 0
                    psiarr(kind,Bind,zind) = 0 ;
                else
                    if Bpr <= Bmin
                        Emuprvec = real(Emuprarr(:,1,zind)) ;
                        % display('Bprime is less than Bmin') ;
                    elseif Bpr >= Bmax
                        Emuprvec = real(Emuprarr(:,Bsize,zind)) ;
                        % display('Bprime is greater than Bmax') ;
                    elseif Bpr > Bmin && Bpr < Bmax
                        [bval,binx] = ...
                            min((Bpr*ones(size(Bvec)) - Bvec).^2) ;
                        blinx = (Bpr >= Bvec(binx))*binx + ...
                                (Bpr < Bvec(binx))*(binx - 1) ;
                        bhinx = (Bpr >= Bvec(binx))*(binx + 1) + ...
                                (Bpr < Bvec(binx))*binx ;
                        blprob = (Bvec(bhinx) - Bpr)/...
                                 (Bvec(bhinx) - Bvec(blinx)) ;
                        bhprob = 1 - blprob ;
                        Emuprvec = real(blprob*Emuprarr(:,blinx,zind) + ...
                                   bhprob*Emuprarr(:,bhinx,zind)) ;
                    end
                    if sum(imag(Emuprvec) ~= 0) > 0
                        psiarr(kind,Bind,zind) = 0 ;
                    else
                        params  = [beta delta gamma kt wt rt taut d] ;
                        kprinit = ((1-taut)*wt + (1+rt-delta)*kt + d)/2 ;
                        [kpr,fval_psi] = fsolve(@unstabEuler_vfi,kprinit,options,...
                                         Emuprvec,kvec,params) ;
                        % [kpr,~] = fmincon(@psiEuler_vfi,kprinit,[],...
                        %    [],[],[],kprmin,kprmax,[],options,params,EVprvec,kvec) ;
                        psiarr(kind,Bind,zind) = kpr ;
                        if fval_psi > 0.1
                            display([psiiter,kind,Bind,zind]) ;
                            error('Euler equation did not solve') ;
                        end
                    end
                end
            end
        end
    end
    psidist = sum((psiprime(:) - psiarr(:)).^2) ;
    disp(psiiter) ;
    disp(psidist) ;
end

%-------------------------------------------------------------------------%
% Solve for value function without optimization using psimat
%-------------------------------------------------------------------------%
% cmat = wmat + ((1 - delta)*ones(size(rmat)) + rmat).*kmat - psimat ;
% if gamma == 1
%     umat = log(cmat) ;
% elseif gamma > 1
%     umat = (cmat.^(1-gamma) - ones(size(cmat)))./(1-gamma) ;
% end
% Vinit = umat ;
% Vnew = Vinit ;
% 
% % Start iteration loop here
% Vpr = Vnew ;
% EVprmat = sum(repmat(reshape(Vpr,[ksize,1,zsize]),[1,zsize,1]).*...
%               repmat(reshape(zTrans,[1,zsize,zsize]),[ksize,1,1]),3) ;
% EVpr = ...
% Vnew = umat + beta*EVpr ;

%-------------------------------------------------------------------------%
% Generate two IRFs using VFI policy function k' = psi(k,z)
%-------------------------------------------------------------------------%
% kirfvec_vfi1 = 1 x (Tirf+1) vector, IRF1 path for capital K using state
%                at t=0 of (k_0,z_0) = (kss,mu)
% kirfvec_vfi2 = 1 x (Tirf+1) vector, IRF2 path for capital K using state
%                at t=0 of (k_0,z_0) = (1.1*kss,mu). Note: we have to force
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
% kminind_1    = scalar, index of kvec value to which ktp1_1 is closest
% kminind_2    = scalar, index of kvec value to which ktp1_2 is closest
% kuw          = scalar, upper bound weight for linear interpolation of
%                psimat for integration
% klw          = scalar, lower bound weight for linear interpolation of
%                psimat for integration
% ktp2v_1      = 1 x znodes vector of potential values for ktp2 next period
%                for IRF1
% ktp2v_2      = 1 x znodes vector of potential values for ktp2 next period
%                for IRF2
% ctp1v_1      = 1 x znodes vector of potential values for ctp1 next period
%                for IRF1
% ctp1v_2      = 1 x znodes vector of potential values for ctp1 next period
%                for IRF2
% Emutp1_1     = scalar, discounted expected marginal utility next period
%                for IRF1
% Emutp1_2     = scalar, discounted expected marginal utility next period
%                for IRF2
%-------------------------------------------------------------------------%
% kirfvec_vfi1 = [kss,zeros(1,Tirf)] ;
% kirfvec_vfi2 = [1.1*kss,zeros(1,Tirf)] ;
% eulvec_vfi1  = zeros(1,Tirf) ;
% eulvec_vfi2  = zeros(1,Tirf) ;
% for tind = 1:Tirf
%     zt_1 = zirfvec1(tind) ;
%     zt_2 = zirfvec2(tind) ;
%     kt_1 = kirfvec_vfi1(tind) ;
%     kt_2 = kirfvec_vfi2(tind) ;
% 
%     if zt_1 <= zvec(1) && kt_1 <= kvec(1)
%         ktp1_1 = psimat(1,1) ;
%     elseif zt_1 >= zvec(zsize) && kt_1 >= kvec(ksize)
%         ktp1_1 = psimat(ksize,zsize) ;
%     elseif zt_1 <= zvec(1) && kt_1 >= kvec(ksize)
%         ktp1_1 = psimat(ksize,1) ;
%     elseif zt_1 >= zvec(zsize) && kt_1 <= kvec(1)
%         ktp1_1 = psimat(1,zsize) ;
%     elseif zt_1 <= zvec(1) && kt_1 > kvec(1) && kt_1 < kvec(ksize)
%         ktp1_1 = interp1(kvec,psimat(:,1),kt_1,'spline') ;
%     elseif zt_1 >= zvec(zsize) && kt_1 > kvec(1) && kt_1 < kvec(ksize)
%         ktp1_1 = interp1(kvec,psimat(:,zsize),kt_1,'spline') ;
%     elseif zt_1 > zvec(1) && zt_1 < zvec(zsize) && kt_1 <= kvec(1)
%         ktp1_1 = interp1(zvec,psimat(1,:),zt_1,'spline') ;
%     elseif zt_1 > zvec(1) && zt_1 < zvec(zsize) && kt_1 >= kvec(ksize)
%         ktp1_1 = interp1(zvec,psimat(ksize,:),zt_1,'spline') ;
%     elseif zt_1 > zvec(1) && zt_1 < zvec(zsize) && ...
%       kt_1 > kvec(1) && kt_1 < kvec(ksize)
%         ktp1_1 = interp2(zmat,kmat,psimat,zt_1,kt_1,'spline') ;
%     end
% 
%     if zt_2 <= zvec(1) && kt_2 <= kvec(1)
%         ktp1_2 = psimat(1,1) ;
%     elseif zt_2 >= zvec(zsize) && kt_2 >= kvec(ksize)
%         ktp1_2 = psimat(ksize,zsize) ;
%     elseif zt_2 <= zvec(1) && kt_2 >= kvec(ksize)
%         ktp1_2 = psimat(ksize,1) ;
%     elseif zt_2 >= zvec(zsize) && kt_2 <= kvec(1)
%         ktp1_2 = psimat(1,zsize) ;
%     elseif zt_2 <= zvec(1) && kt_2 > kvec(1) && kt_2 < kvec(ksize)
%         ktp1_2 = interp1(kvec,psimat(:,1),kt_2,'spline') ;
%     elseif zt_2 >= zvec(zsize) && kt_2 > kvec(1) && kt_2 < kvec(ksize)
%         ktp1_2 = interp1(kvec,psimat(:,zsize),kt_2,'spline') ;
%     elseif zt_2 > zvec(1) && zt_2 < zvec(zsize) && kt_2 <= kvec(1)
%         ktp1_2 = interp1(zvec,psimat(1,:),zt_2,'spline') ;
%     elseif zt_2 > zvec(1) && zt_2 < zvec(zsize) && kt_2 >= kvec(ksize)
%         ktp1_2 = interp1(zvec,psimat(ksize,:),zt_2,'spline') ;
%     elseif zt_2 > zvec(1) && zt_2 < zvec(zsize) && ...
%       kt_2 > kvec(1) && kt_2 < kvec(ksize)
%         ktp1_2 = interp2(zmat,kmat,psimat,zt_2,kt_2,'spline') ;
%     end
% 
%     kirfvec_vfi1(tind+1) = ktp1_1 ;
%     kirfvec_vfi2(tind+1) = ktp1_2 ;
%     wt_1 = (1-alpha)*exp(zt_1)*(kt_1^alpha) ;
%     wt_2 = (1-alpha)*exp(zt_2)*(kt_2^alpha) ;
%     rt_1 = alpha*exp(zt_1)*(kt_1^(alpha-1)) ;
%     rt_2 = alpha*exp(zt_2)*(kt_2^(alpha-1)) ;
%     ct_1 = wt_1 + (1 + rt_1 - delta)*kt_1 - ktp1_1 ;
%     ct_2 = wt_2 + (1 + rt_2 - delta)*kt_2 - ktp1_2 ;
%     mut_1 = ct_1^(-gamma) ;
%     mut_2 = ct_2^(-gamma) ;
% 
%     [~,zminind_1] = min((repmat(zt_1,[length(zvec),1]) - zvec).^2,[],1) ;
%     if zt_1 <= zvec(1)
%         zprobs_1 = zTrans(1,:) ;
%     elseif zt_1 <= zvec(zminind_1)
%         zuw_1    = (zt_1 - zvec(zminind_1-1))/...
%                    (zvec(zminind_1) - zvec(zminind_1-1)) ;
%         zlw_1    = 1 - zuw_1 ;
%         zprobs_1 = zlw_1*zTrans(zminind_1-1,:) + ...
%                    zuw_1*zTrans(zminind_1,:) ;
%     elseif zt_1 > zvec(zminind_1) && zt_1 < zvec(length(zvec))
%         zlw_1    = (zvec(zminind_1+1) - zt_1)/(zvec(zminind_1+1) - ...
%                    zvec(zminind_1)) ;
%         zuw_1    = 1 - zlw_1 ;
%         zprobs_1 = zlw_1*zTrans(zminind_1,:) + ...
%                    zuw_1*zTrans(zminind_1+1,:) ;
%     elseif zt_1 >= zvec(length(zvec))
%         zprobs_1 = zTrans(length(zvec),:) ;
%     end
%     rtp1v_1 = (alpha*ktp1_1^(alpha-1)).*exp(zvec') ;
%     wtp1v_1 = ((1-alpha)*(ktp1_1^alpha)).*exp(zvec') ;
%     [~,kminind_1] = min(abs(ktp1_1 - kvec)) ;
%     if ktp1_1 <= kvec(1)
%         ktp2v_1 = psimat(1,:) ;
%     elseif ktp1_1 >= kvec(ksize)
%         ktp2v_1 = psimat(ksize,:) ;
%     elseif ktp1_1 > kvec(1) && ktp1_1 < kvec(kminind_1)
%         kuw = (ktp1_1 - kvec(kminind_1-1))/...
%               (kvec(kminind_1) - kvec(kminind_1-1)) ;
%         klw = 1 - kuw ;
%         ktp2v_1 = kuw*psimat(kminind_1,:) + klw*psimat(kminind_1-1,:) ;
%     elseif ktp1_1 < kvec(ksize) && ktp1_1 > kvec(kminind_1)
%         kuw = (ktp1_1 - kvec(kminind_1))/...
%               (kvec(kminind_1+1) - kvec(kminind_1)) ;
%         klw = 1 - kuw ;
%         ktp2v_1 = kuw*psimat(kminind_1+1,:) + klw*psimat(kminind_1,:) ;
%     end
%     ctp1v_1 = wtp1v_1 + ktp1_1*((1-delta)*ones(1,zsize) + rtp1v_1) - ...
%               ktp2v_1 ;
%     Emutp1_1 = sum(zprobs_1.*((1-delta)*ones(1,zsize) + rtp1v_1).*...
%                (ctp1v_1.^(-gamma))) ;
%     eulvec_vfi1(tind) = beta*Emutp1_1/mut_1 - 1 ;
% 
%     [~,zminind_2] = min((repmat(zt_2,[length(zvec),1]) - zvec).^2,[],1) ;
%     if zt_2 <= zvec(1)
%         zprobs_2 = zTrans(1,:) ;
%     elseif zt_2 <= zvec(zminind_2)
%         zuw_2    = (zt_2 - zvec(zminind_2-1))/...
%                    (zvec(zminind_2) - zvec(zminind_2-1)) ;
%         zlw_2    = 1 - zuw_2 ;
%         zprobs_2 = zlw_2*zTrans(zminind_2-1,:) + ...
%                    zuw_2*zTrans(zminind_2,:) ;
%     elseif zt_2 > zvec(zminind_2) && zt_2 < zvec(length(zvec))
%         zlw_2    = (zvec(zminind_2+1) - zt_2)/(zvec(zminind_2+1) - ...
%                    zvec(zminind_2)) ;
%         zuw_2    = 1 - zlw_2 ;
%         zprobs_2 = zlw_2*zTrans(zminind_2,:) + ...
%                    zuw_2*zTrans(zminind_2+1,:) ;
%     elseif zt_2 >= zvec(length(zvec))
%         zprobs_2 = zTrans(length(zvec),:) ;
%     end
%     rtp1v_2 = (alpha*ktp1_2^(alpha-1)).*exp(zvec') ;
%     wtp1v_2 = ((1-alpha)*(ktp1_2^alpha)).*exp(zvec') ;
%     [~,kminind_2] = min(abs(ktp1_2 - kvec)) ;
%     if ktp1_2 <= kvec(1)
%         ktp2v_2 = psimat(1,:) ;
%     elseif ktp1_2 >= kvec(ksize)
%         ktp2v_2 = psimat(ksize,:) ;
%     elseif ktp1_2 > kvec(1) && ktp1_2 < kvec(kminind_2)
%         kuw = (ktp1_2 - kvec(kminind_2-1))/...
%               (kvec(kminind_2) - kvec(kminind_2-1)) ;
%         klw = 1 - kuw ;
%         ktp2v_2 = kuw*psimat(kminind_2,:) + klw*psimat(kminind_2-1,:) ;
%     elseif ktp1_2 < kvec(ksize) && ktp1_2 > kvec(kminind_2)
%         kuw = (ktp1_2 - kvec(kminind_2))/...
%               (kvec(kminind_2+1) - kvec(kminind_2)) ;
%         klw = 1 - kuw ;
%         ktp2v_2 = kuw*psimat(kminind_2+1,:) + klw*psimat(kminind_2,:) ;
%     end
%     ctp1v_2 = wtp1v_2 + ktp1_2*((1-delta)*ones(1,zsize) + rtp1v_2) - ...
%               ktp2v_2 ;
%     Emutp1_2 = sum(zprobs_2.*((1-delta)*ones(1,zsize) + rtp1v_2).*...
%                (ctp1v_2.^(-gamma))) ;
%     eulvec_vfi2(tind) = beta*Emutp1_2/mut_2 - 1 ;
% end
% 
% %-------------------------------------------------------------------------%
% % Generate simulated time path for K using VFI policy function
% % k' = psi(k,z)
% %-------------------------------------------------------------------------%
% % ksimvec_vfi3 = 1 x Tsim+1 vector, simulated path of capital K starting at
% %                state (k_0,z_0) = (kss,mu)
% % eulvec_vfi3  = 1 x Tsim vector of Euler errors from longer simulation
% % zt_3         = scalar, current period z value for simulation
% % kt_3         = scalar, current period k value for simulation
% % ktp1_3       = scalar, next period capital stock value for simulation
% % wt_3         = scalar, current period w value for simulation
% % rt_3         = scalar, current period r value for simulation
% % ct_3         = scalar, current period c value for simulation
% % mut_3        = scalar, current period marginal utility for simulation
% % zminind_3    = integer, index of value in discretized support of z to
% %                which zt_3 is closest
% % zprobs_3     = 1 x znodes vector, probability weights for ztp1_3 tomorrow
% %                which are a linear interpolation of the transition matrix
% %                row for which zt_3 is closest
% % zuw_3        = scalar, weight on upper node of zvec for convex
% %                combination of transition matrix rows for simulation
% % zlw_3        = scalar, weight on lower node of zvec for convex
% %                combination of transition matrix rows for simulation
% % rtp1v_3      = 1 x znodes vector of potential values for rtp1 next period
% %                for simulation
% % wtp1v_3      = 1 x znodes vector of potential values for wtp1 next period
% %                for simulation
% % kminind_3    = scalar, index of kvec value to which ktp1_3 is closest
% % ktp2v_3      = 1 x znodes vector of potential values for ktp2 next period
% %                for simulation
% % ctp1v_3      = 1 x znodes vector of potential values for ctp1 next period
% %                for simulation
% % Emutp1_3     = scalar, discounted expected marginal utility next period
% %                for simulation
% %-------------------------------------------------------------------------%
% % Fit a polynomial k'=exp(a0+a1*k+a2*z+a3*k^2+a4*z^2+a5*k*z) to psimat in
% % order to interpolate
% ktp1vec = reshape(psimat,[ksize*zsize,1]) ;
% ktvec = reshape(repmat(kvec,[1,zsize]),[ksize*zsize,1]) ;
% ztvec = reshape(repmat(zvec',[ksize,1]),[ksize*zsize,1]) ;
% xx = [ones(size(ktvec)),ktvec,ztvec,ktvec.^2,ztvec.^2,ktvec.*ztvec] ;
% psicoefs = (xx'*xx)\(xx'*ktp1vec) ;
% 
% ksimvec_vfi3 = [kss,zeros(1,Tsim)] ;
% eulvec_vfi3  = zeros(1,Tsim) ;
% for tind = 1:Tsim
%     zt_3 = zsimvec3(tind) ;
%     kt_3 = ksimvec_vfi3(tind) ;
%     if zt_3 <= zvec(1) && kt_3 <= kvec(1)
%         ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
%     elseif zt_3 >= zvec(zsize) && kt_3 >= kvec(ksize)
%         ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
%     elseif zt_3 <= zvec(1) && kt_3 >= kvec(ksize)
%         ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
%     elseif zt_3 >= zvec(zsize) && kt_3 <= kvec(1)
%         ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
%     elseif zt_3 <= zvec(1) && kt_3 > kvec(1) && kt_3 < kvec(ksize)
%         ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
%     elseif zt_3 >= zvec(zsize) && kt_3 > kvec(1) && kt_3 < kvec(ksize)
%         ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
%     elseif zt_3 > zvec(1) && zt_3 < zvec(zsize) && kt_3 <= kvec(1)
%         ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
%     elseif zt_3 > zvec(1) && zt_3 < zvec(zsize) && kt_3 >= kvec(ksize)
%         ktp1_3 = max([1,kt_3,zt_3,kt_3^2,zt_3^2,kt_3*zt_3]*psicoefs,10^(-6)) ;
%     elseif zt_3 > zvec(1) && zt_3 < zvec(zsize) && ...
%       kt_3 > kvec(1) && kt_3 < kvec(ksize)
%         ktp1_3 = interp2(zmat,kmat,psimat,zt_3,kt_3,'spline') ;
%     end
%     ksimvec_vfi3(tind+1) = ktp1_3 ;
%     wt_3 = (1-alpha)*exp(zt_3)*(kt_3^alpha) ;
%     rt_3 = alpha*exp(zt_3)*(kt_3^(alpha-1)) ;
%     ct_3 = wt_3 + (1 + rt_3 - delta)*kt_3 - ktp1_3 ;
%     mut_3 = ct_3^(-gamma) ;
% 
%     [~,zminind_3] = min((repmat(zt_3,[length(zvec),1]) - zvec).^2,[],1) ;
%     if zt_3 <= zvec(1)
%         zprobs_3 = zTrans(1,:) ;
%     elseif zt_3 >= zvec(length(zvec))
%         zprobs_3 = zTrans(length(zvec),:) ;
%     elseif zt_3 <= zvec(zminind_3)
%         zuw_3    = (zt_3 - zvec(zminind_3-1))/...
%                    (zvec(zminind_3) - zvec(zminind_3-1)) ;
%         zlw_3    = 1 - zuw_3 ;
%         zprobs_3 = zlw_3*zTrans(zminind_3-1,:) + ...
%                    zuw_3*zTrans(zminind_3,:) ;
%     elseif zt_3 > zvec(zminind_3) && zt_3 < zvec(length(zvec))
%         zlw_3    = (zvec(zminind_3+1) - zt_3)/(zvec(zminind_3+1) - ...
%                    zvec(zminind_3)) ;
%         zuw_3    = 1 - zlw_3 ;
%         zprobs_3 = zlw_3*zTrans(zminind_3,:) + ...
%                    zuw_3*zTrans(zminind_3+1,:) ;
%     end
%     rtp1v_3 = (alpha*ktp1_3^(alpha-1)).*exp(zvec') ;
%     wtp1v_3 = ((1-alpha)*(ktp1_3^alpha)).*exp(zvec') ;
% %     ktp2v_3 = max(psicoefs(1)*ones(1,zsize) + psicoefs(2)*ktp1_3*ones(1,zsize) ...
% %               + psicoefs(3)*zvec' + (psicoefs(4)*ktp1_3^2)*ones(1,zsize)...
% %               + psicoefs(5)*zvec'.^2 + (psicoefs(6)*ktp1_2)*zvec',10^(-6)) ;
%     [~,kminind_3] = min(abs(ktp1_3 - kvec)) ;
%     if ktp1_3 <= kvec(1)
%         ktp2v_3 = psimat(1,:) ;
%     elseif ktp1_3 >= kvec(ksize)
%         ktp2v_3 = psimat(ksize,:) ;
%     elseif ktp1_3 > kvec(1) && ktp1_3 < kvec(kminind_3)
%         kuw = (ktp1_3 - kvec(kminind_3-1))/...
%               (kvec(kminind_3) - kvec(kminind_3-1)) ;
%         klw = 1 - kuw ;
%         ktp2v_3 = kuw*psimat(kminind_3,:) + klw*psimat(kminind_3-1,:) ;
%     elseif ktp1_3 < kvec(ksize) && ktp1_3 > kvec(kminind_3)
%         kuw = (ktp1_3 - kvec(kminind_3))/...
%               (kvec(kminind_3+1) - kvec(kminind_3)) ;
%         klw = 1 - kuw ;
%         ktp2v_3 = kuw*psimat(kminind_3+1,:) + klw*psimat(kminind_3,:) ;
%     end
%     ctp1v_3 = wtp1v_3 + ktp1_3*((1-delta)*ones(1,zsize) + rtp1v_3) - ...
%               ktp2v_3 ;
%     Emutp1_3 = sum(zprobs_3.*((1-delta)*ones(1,zsize) + rtp1v_3).*...
%                (ctp1v_3.^(-gamma))) ;
%     eulvec_vfi3(tind) = beta*Emutp1_3/mut_3 - 1 ;
% end
% 
% %-------------------------------------------------------------------------%
% % Calculate two Euler error distance measures for the two IRFs and the one
% % simulation
% %-------------------------------------------------------------------------%
% % eulerr_rms_vfi1 = scalar, root mean squared Euler error from IRF1
% % eulerr_mab_vfi1 = scalar, maximum absolute Euler error from IRF1
% % eulerr_rms_vfi2 = scalar, root mean squared Euler error from IRF2
% % eulerr_mab_vfi2 = scalar, maximum absolute Euler error from IRF1
% % eulerr_rms_vfi3 = scalar, root mean squared Euler error from simulation
% % eulerr_mab_vfi3 = scalar, maximum absolute Euler error from simulation
% %-------------------------------------------------------------------------%
% eulerr_rms_vfi1 = sqrt(sum(eulvec_vfi1.^2)/Tirf) ;
% eulerr_mab_vfi1 = max(abs(eulvec_vfi1)) ;
% eulerr_rms_vfi2 = sqrt(sum(eulvec_vfi2.^2)/Tirf) ;
% eulerr_mab_vfi2 = max(abs(eulvec_vfi2)) ;
% eulerr_rms_vfi3 = sqrt(sum(eulvec_vfi3.^2)/Tsim) ;
% eulerr_mab_vfi3 = max(abs(eulvec_vfi3)) ;
% 
% display([eulerr_rms_vfi1,eulerr_mab_vfi1;...
%     eulerr_rms_vfi2,eulerr_mab_vfi2;...
%     eulerr_rms_vfi3,eulerr_mab_vfi3]) ;
% 
% %-------------------------------------------------------------------------%
% % Plot IRFs and simulation
% %-------------------------------------------------------------------------%
% figure(2)
% plot((0:1:Tirf),kirfvec_vfi1,(0:1:Tirf),kirfvec_vfi2) ;
% xlabel('Period (t)','FontSize', 12) ;
% ylabel('Aggregate Capital (K)','FontSize', 12) ;
% 
% figure(3)
% plot((0:1:Tsim),ksimvec_vfi3) ;
% xlabel('Period (t)','FontSize', 12) ;
% ylabel('Aggregate Capital (K)','FontSize', 12) ;

%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save unstab_vfi.mat ;
