% This file is written in order to calculate run value function
% iteration on a standard infinte horizon model.  We will build a grid
% that goes from -3 times the certainty equivalence ss and goes to 3
% times certainty equivalence ss.

% make params global
global cbeta ;
global cgamma ;
global calpha ;
global cdelta ;
global L ;
global cmu ;
global csigma ;

% make necessary vectors/matrices/variables global
global kvec ;
global zvec ;
global kmin ;
global kmax ;
global zmin ;
global zmax ;
global V_init ;
global kt ;
global zt ;

% Here we initialize our parameters
cbeta       = 0.96 ;
cgamma      = 2 ;
calpha      = 0.35 ;
cdelta      = 0.05 ;
% L           = 1 ;
% cmu         = 0 ;
crho_qtr    = 0.95 ;
crho        = crho_qtr^(4) ;
csigma_qtr  = 0.007 ;
crhosigsum  = 0 ;
for sind = 1:4
    crhosigsum = crhosigsum + crho_qtr^(2*(sind-1)) ;
end
csigma   = sqrt(crhosigsum)*csigma_qtr ;

% Solve for the certainty equivalence steady state
rcss = 1/cbeta - 1 + cdelta ;
kcss = (rcss/calpha)^(1/(calpha-1)) ;

% We need to choose the size of each vector and set the ranges for the
% grid that we will be using.  Also, we choose tolerance and max iters
% for convergence here.
zsize   = 15 ;
ksize   = 300 ;
kmin    = 1e-5 ;
kmax    = 3 * kcss ;
zmin    = -5 * csigma ;
zmax    = 5 * csigma ;
kvec    = linspace(kmin, kmax, ksize) ;
zvec    = linspace(zmin, zmax, zsize) ;
tol     = 1e-8 ;
maxiter = 1500 ;
optionsp = optimset('Display','off','MaxFunEvals',100000,'MaxIter',1000,...
                    'TolFun',1e-15,'Algorithm','active-set') ;

% Set an initial value function and initialize several other variables
% necessary for the VFI
V_init   = zeros(ksize, zsize) ;
V_prime  = zeros(ksize, zsize) ;
kpchoice = zeros(ksize, zsize) ;
V_diff   = 10 ;
iters    = 0 ;


while V_diff > tol & iters < maxiter ;
    tic
    iters = iters + 1 ;
    V_init = V_prime ;
    disp(iters) ;
    disp(V_diff) ;
    for kind = 1:ksize ;
        for zind = 1:zsize ;
            kt = kvec(kind) ;
            zt = zvec(zind) ;
            ktp1 = fmincon(@psisolve2,kt,[],[],[],...
               [],10^(-10),kmax - 2, [],optionsp) ;

            % ktp1 = fmincon(@psisolve2, k, kmin, kmax) ;
            V_prime(kind, zind) = -psisolve2(ktp1) ;
            kpchoice(kind, zind) = ktp1 ;
        end
    end
    V_diff = norm(V_prime - V_init, 2) ;
    toc
end




