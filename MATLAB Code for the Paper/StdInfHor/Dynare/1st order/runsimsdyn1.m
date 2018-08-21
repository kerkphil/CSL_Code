% Dynare Eulers 1st order
% We want to run Dynare and then bring out the policy functions and do 
% numerical integration to calculate the true Euler Errors from our model

% Define number of periods we will simulate over

% First we run our model in Dynare
dynare InfHorDyn.mod noclearall
T = 150
sigma = 0.0130058759370760;

% Display order of variable declartion so we can define policy funcs
disp('Note the importance of order of variable declaration')
disp('w, Y, K, z, c, r')

% Gives steady states b/c First order approx in Dynare is done by
% y_t = y* + A y_(t-1) + B u_t where y* is ss, yh_(t-1) is the set of values
% from the period before - the ss (y(t-1) - y*)
% period before and u_t are the shocks from this period.  A, B are the
% policy function values.
ystar = oo_.dr.ys ;
A = oo_.dr.ghx ;
B = oo_.dr.ghu ;
mu = 0 ;

% Gives us the matrix we need to simulate for 150 periods
load('simulationshocks.mat');
save('zshockdynirf1', 'z_eps')
save('params', 'cbeta', 'cgamma', 'cdelta', 'mu', 'csigma', 'crho', 'calpha', 'ystar', 'A', 'B')

% Initializes matrix of variables and ss values
varmatrix = zeros(6,T) ;
varmatrix(:, 1) = ystar ;
kmatsim = zeros(T, 1) ;
zmatsim = zeros(T, 1) ;
kss = ystar(3) ;
zss = ystar(4) ;
kmatsim(1, 1) = ystar(3) ;
zmatsim(1, 1) = ystar(4) ;

%  Stochastic Simulation
for i = 2:T
    % Get the values of k and z from last period to plug into policy funcs
    kztemp = [kmatsim(i - 1, 1) - kss; zmatsim(i - 1, 1) - zss] ;
    
    % Here we plug in the values necessary to evaluate the policy function
    % and find the values of our variables for t+1
    varmatrix(:, i) = ystar + A * kztemp + B .* simeps(i) ;
    kmatsim(i, 1) = varmatrix(3, i) ;
    zmatsim(i, 1) = varmatrix(4, i) ;
end

% Save vectors so we can plot the K series and compare
save('kmatsimseries1ord.mat', 'kmatsim', 'zmatsim', 'simeps', 'kss')

save('kmatirf1ord.mat', 'K_eps', 'kss', 'z_eps')

% This section contains the code to give the model a 2
% standard deviation shock to productivity.
% Initializes matrix of variables and ss values
var2matrix = zeros(6,T) ;
var2matrix(:, 1) = ystar ;
var2matrix(3, 1) = 2.0 * kss ;
kmat2 = zeros(T, 1) ;
zmat2 = zeros(T ,1) ;
kmat2(1, 1) = ystar(3)* 2.0 ;
zmat2(1, 1) = ystar(4) ;
shock2eps2 = zeros(T, 1) ;
shock2eps2(2, 1) = 2.0 * sigma ;

for i = 2:T
    % Get the values of k and z from last period to plug into policy funcs
    kztemp2 = [kmat2(i - 1, 1) - kss; zmat2(i - 1, 1) - zss] ;
    
    % Here we plug in the values necessary to evaluate the policy function
    % and find the values of our variables for t+1
    var2matrix(:, i) = ystar + A * kztemp2 + B .* shock2eps2(i) ;
    kmat2(i, 1) = var2matrix(3, i) ;
    zmat2(i, 1) = var2matrix(4, i) ;
end

save('kmat2irf1ord.mat', 'kmat2', 'zmat2', 'shock2eps2', 'kss')