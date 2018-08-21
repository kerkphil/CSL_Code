function diff = vMtaxEuler_vfi(kpr,Emuprvec,kvec,params)
%--------------------------------------------------------------------------
% kpr      = scalar, next period capital stock
% Emuprvec = ksize x 1 vector, the column of Emuprmat E[mu(k',z)]
%            corresponding to zind
% kvec     = ksize x 1 vector of grid points in k dimension
% params   = 1 x 10 vector, parameters passed in to program
% beta     = scalar, discount factor
% delta    = scalar, depreciation rate
% gamma    = scalar, coefficient of relative risk aversion
% kt       = scalar, current period capital stock
% wt       = scalar, current period real wage
% rt       = scalar, current period real interest (rental) rate
% tau1     = scalar, asymptotic marginal tax rate as income goes to to
%            negative infinity
% tau2     = scalar, asymptotic marginal tax rate as income goes to to
%            positive infinity
% atau     = scalar, transition parameter in marginal tax rate function
% btau     = scalar, left/right shifter in marginal tax rate function
%--------------------------------------------------------------------------
beta = params(1) ;
delta = params(2) ;
gamma = params(3) ;
kt = params(4) ;
wt = params(5) ;
rt = params(6) ;
tau1 = params(7) ;
tau2 = params(8) ;
atau = params(9) ;
btau = params(10) ;

%--------------------------------------------------------------------------
% Calculate the error (difference) version of the Euler equation
%--------------------------------------------------------------------------
% ct  = scalar, current period consumption
% mut = scalar, current period marginal utility of consumption
% m   = scalar, slope parameter for linear interpolation
% b   = scalar, intercept parameter for linear interpolation
%--------------------------------------------------------------------------
xt = wt + (rt - delta)*kt ;
taut = tau1 + ((1/pi)*atan(atau*xt + btau) + 0.5)*(tau2 - tau1) ;
dt = taut*xt ;
ct = wt + (1 + rt - delta)*kt - taut*xt - kpr + dt ;
mut = ct^(-gamma) ;

% Approximate value for Emupr with cubic spline interpolation on interior
% of Emuprvec and linear extrapolation on exterior of Emuprvec
if kpr <= kvec(1)
    % linear extrapolation from lowest two points
    m = (Emuprvec(2) - Emuprvec(1))/(kvec(2) - kvec(1)) ;
    b = -m*kvec(1) + Emuprvec(1) ;
    Emupr = m*kpr + b ;
elseif kpr > kvec(1) && kpr < kvec(length(kvec))
    % cubic spline interpolation between points on interior of grid
    Emupr = interp1(kvec,Emuprvec,kpr,'spline') ;
elseif kpr >= kvec(length(kvec)) ;
    % linear extrapolation from highest two points
    m = (Emuprvec(length(Emuprvec)) - Emuprvec(length(Emuprvec)-1))/...
        (kvec(length(kvec)) - kvec(length(kvec)-1)) ;
    b = -m*kvec(length(kvec)) + Emuprvec(length(Emuprvec)) ;
    Emupr = m*kpr + b ;
end

diff = mut - beta*Emupr ;
