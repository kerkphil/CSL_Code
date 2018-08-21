function diff = unstabEuler_vfi(kpr,Emuprvec,kvec,params)
%--------------------------------------------------------------------------
% kpr      = scalar, next period capital stock
% Emuprvec = ksize x 1 vector, the column of Emuprmat E[mu(k',B',z)]
%            corresponding to B' and zind
% kvec     = ksize x 1 vector of grid points in k dimension
% params   = 1 x 8 vector, parameters passed in to program
% beta     = scalar, discount factor
% delta    = scalar, depreciation rate
% gamma    = scalar, coefficient of relative risk aversion
% kt       = scalar, current period capital stock
% wt       = scalar, current period real wage
% rt       = scalar, current period real interest (rental) rate
% Bt       = scalar, current period value of B
% taut     = scalar, current period value of tau
% d        = scalar, fixed lump-sum transfer
%--------------------------------------------------------------------------
beta = params(1) ;
delta = params(2) ;
gamma = params(3) ;
kt = params(4) ;
wt = params(5) ;
rt = params(6) ;
taut = params(7) ;
d = params(8) ;

%--------------------------------------------------------------------------
% Calculate the error (difference) version of the Euler equation
%--------------------------------------------------------------------------
% ct  = scalar, current period consumption
% mut = scalar, current period marginal utility of consumption
% m   = scalar, slope parameter for linear interpolation
% b   = scalar, intercept parameter for linear interpolation
%--------------------------------------------------------------------------
ct = (1 - taut)*wt + (1 + rt - delta)*kt - kpr + d ;
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
