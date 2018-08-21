% Dynare Eulers 2nd order
% We want to run Dynare and then bring out the policy functions and do 
% numerical integration to calculate the true Euler Errors from our model

% Define number of periods we will simulate over
% Small Change for spencer
% First we run our model in Dynare
dynare InfHorDyn2.mod noclearall
T = 150 ;
sigma = 0.0130058759370760 ;

%Display order of variable declartion so we can define policy funcs
disp('Note the importance of order of variable declaration')
disp('w, Y, K, z, c, r')

% Gives steady states b/c Second order approx in Dynare is done by
% y_t = y* + delta^2 + A y_t-1 + B u_t + .5 C (y_t-1 kr y_t-1) + .5 D (u_t kr u_t) +
% E (y_t-1 kr u_t) where y* is ss, yh_(t-1) is the set of values from the 
% period before - the ss (y(t-1) - y*) period before and u_t are the shocks
% from this period.  A, B, C, D, E, deltasqr are the policy function values.
ystar = oo_.dr.ys ;
deltasqr = oo_.dr.ghs2 ;
A = oo_.dr.ghx ;
B = oo_.dr.ghu ;
C = oo_.dr.ghxx ;
D = oo_.dr.ghuu ;
E = oo_.dr.ghxu ;
mu = 0 ;
kss2 = ystar(3) + deltasqr(3) ; % Maybe add the dltsqr term for k in here?

save('params', 'cbeta', 'cgamma', 'cdelta', 'mu', 'csigma', 'crho', ...
    'calpha', 'ystar', 'deltasqr', 'A', 'B', 'C', 'D', 'E')
% Gives us the matrix we need to simulate for 150 periods
load('simulationshocks.mat');
save('kmatirf2ord.mat', 'K_eps', 'kss2')

% Initializes matrix of variables and ss values
varmatrix = zeros(6,T) ;
varmatrix(:, 1) = ystar ;
kmatsim = zeros(T, 1) ;
zmatsim = zeros(T ,1) ;
zss = ystar(4) ;
kmatsim(1, 1) = ystar(3) ;
zmatsim(1, 1) = ystar(4) ;

for i = 2:T
    %Get the values of k and z from last period to plug into policy funcs
    kk_eps = kmatsim(i - 1, 1) - kss2 ;
    ztm_eps = zmatsim(i - 1, 1) - zss ;
    kztempkron = [kk_eps^2, kk_eps*ztm_eps, kk_eps*ztm_eps, ztm_eps^2]' ;
    kztemp = [kmatsim(i - 1, 1) - kss2; zmatsim(i - 1, 1) - zss] ;
    eps_tt = simeps(i) ;
    
    
    %Here we plug in the values necessary to evaluate the policy function
    %and find the values of our variables for t+1
    varmatrix(:, i) = (ystar + .5 * deltasqr + A*kztemp + .5 .* (C*kztempkron))...
    + B .* eps_tt + .5 .* D .* eps_tt.^2 + E * kztemp .* eps_tt ;
    kmatsim(i, 1) = varmatrix(3, i) ;
    zmatsim(i, 1) = varmatrix(4, i) ;
end

%Save vectors so we can plot the K series and compare
save('kmatsimseries2ord.mat', 'kmatsim', 'zmatsim', 'simeps')


% This section contains the code to give the model a 2
% standard deviation shock to productivity.
% Initializes matrix of variables and ss values
var2matrix = zeros(6,T) ;
var2matrix(:, 1) = ystar ;
kmat22 = zeros(T, 1) ;
zmat22 = zeros(T ,1) ;
kmat22(1, 1) = ystar(3)* 2.0 ;
var2matrix(3, 1) = ystar(3) * 2.0 ;
zmat22(1, 1) = ystar(4) ;
shock2eps = zeros(T, 1) ;
shock2eps(2) = 2.0 * sigma ;

for i = 2:T
    %Get the values of k and z from last period to plug into policy funcs
    kk_eps = kmat22(i - 1, 1) - kss2 ;
    ztm_eps = zmat22(i - 1, 1) - zss ;
    kztempkron = [kk_eps^2, kk_eps*ztm_eps, kk_eps*ztm_eps, ztm_eps^2]' ;
    kztemp = [kmat22(i - 1, 1) - kss2; zmat22(i - 1, 1) - zss] ;
    eps_tt = shock2eps(i) ;
    
    
    %Here we plug in the values necessary to evaluate the policy function
    %and find the values of our variables for t+1
    var2matrix(:, i) = (ystar + .5 .* deltasqr + A*kztemp + .5 .* (C*kztempkron))...
    + B.*eps_tt + .5 .* D .* eps_tt.^2 + E * kztemp .* eps_tt ;
    kmat22(i, 1) = var2matrix(3, i) ;
    zmat22(i, 1) = var2matrix(4, i) ;
end

save('kmat2irf2ord.mat', 'kmat22', 'zmat22', 'shock2eps')
