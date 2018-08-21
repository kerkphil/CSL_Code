function val = psisolve2(k)

% This function evaluates the value of the value function at time t

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

% Anonymous function to calculate the value of utility
util = @(c) ((c^(1 - cgamma) - 1) / (1-cgamma)) ;

% We will use linear interpolation so we need to know the indices above
% and below k.
kindlow  = max(sum(k>kvec),1) ;
kindhigh = kindlow + 1 ;

% we use a linear equation (low is gridpoint below k and high is gridpoint
% above k.  Then v(k, :) = v0(low, :) + (k - k(low)) *
%                              (V0(high, :) - V0(low, :))/(k(high) - k(low))
val = V_init(kindlow,:) + (k - kvec(kindlow)) *...
      (V_init(kindhigh,:) - V_init(kindlow,:)) / (kvec(kindhigh) - kvec(kindlow));

% val = interp1(kvec, V_init, k, 'linear') ;

% Solve for consumption to make sure that k makes sense.  If consumption
% is negative then something is wrong so send val to a large neg number
wt = (1 - calpha) * exp(zt) * kt^calpha ;
rt = (calpha) * exp(zt) * kt^(calpha - 1) ;
ct = wt + (1 + rt - cdelta) * kt - k ;

if ct <= 0
    val = -9999999 - 999 * abs(ct) ;
else
    val = util(ct) + cbeta * (val * pdf('norm', zvec, 0, csigma)') ;

end


% We are running a minimization function so we want to minimize the
% negative value of val (aka maximize the original value of val)
val = -abs(val) ;

return