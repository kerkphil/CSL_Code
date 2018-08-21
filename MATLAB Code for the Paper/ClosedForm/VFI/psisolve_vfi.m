function minfunc = psisolve_vfi(kpr,parameters,EVprvec,kvec)
%-------------------------------------------------------------------------%
% Solve for kpr that maximizes the Bellman equation for particular (k,z)
%-------------------------------------------------------------------------%
% kpr        = scalar, value of k'
% parameters = 1 x 4 vector of parameters passed in
% EVprvec    = ksize x 1 vector, values of E[V(k',z)] for all k' and
%              particular z
% kvec       = ksize x 1 vector, discretized support of k or k'
% beta       = scalar, discount factor
% kt         = scalar, current period capital
% zt         = scalar, current period productivity
% wt         = scalar, current period real wage
% rt         = scalar, current period interest rate
% ct         = scalar, current period consumption
% uct        = scalar, utility of current period consumption
% m          = scalar, slope coefficient of extrapolation line
% b          = scalar, y-intercept of extrapolation line
% EVpr       = scalar, interpolated value of E[V(k',z)]
% Vnew       = scalar, value of being at current state V(k,z) for particular
%              choice of k'
%-------------------------------------------------------------------------%
beta  = parameters(1) ;
kt    = parameters(2) ;
wt    = parameters(3) ;
rt    = parameters(4) ;

ct = wt + rt*kt - kpr ;
uct = log(ct) ;

% Approximate grid for value function with linear interpolation
if kpr <= kvec(1)
    % linear extrapolation from lowest two points
    m = (EVprvec(2) - EVprvec(1))/(kvec(2) - kvec(1)) ;
    b = -m*kvec(1) + EVprvec(1) ;
    EVpr = m*kpr + b ;
elseif kpr > kvec(1) && kpr < kvec(length(kvec))
    % cubic spline interpolation between points on interior of grid
    EVpr = interp1(kvec,EVprvec,kpr,'spline') ;
elseif kpr >= kvec(length(kvec)) ;
    % linear extrapolation from highest two points
    m = (EVprvec(length(EVprvec)) - EVprvec(length(EVprvec)-1))/...
        (kvec(length(kvec)) - kvec(length(kvec)-1)) ;
    b = -m*kvec(length(kvec)) + EVprvec(length(EVprvec)) ;
    EVpr = m*kpr + b ;
end

Vnew    = uct + beta*EVpr ;
minfunc = -Vnew ;