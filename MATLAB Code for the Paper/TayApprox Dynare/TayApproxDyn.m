%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
clear global
tic;
global M_ oo_ options_ ys0_ ex0_
options_ = [];
M_.fname = 'TayApproxDyn';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'TayApproxDyn.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'eps';
M_.exo_names_tex = 'eps';
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.param_names = 'alpha';
M_.param_names_tex = 'alpha';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'lbar');
M_.param_names_tex = char(M_.param_names_tex, 'lbar');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names = char(M_.param_names, 'mu');
M_.param_names_tex = char(M_.param_names_tex, 'mu');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 7;
M_.param_nbr = 8;
M_.orig_endo_nbr = 7;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('TayApproxDyn_dynamic');
M_.lead_lag_incidence = [
 0 3 0;
 0 4 10;
 1 5 0;
 0 6 0;
 0 7 11;
 0 8 0;
 2 9 0;]';
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(8, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 23;
M_.NNZDerivatives(2) = 18;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.3;
alpha = M_.params( 1 );
M_.params( 2 ) = 0.99;
beta = M_.params( 2 );
M_.params( 3 ) = 0.025;
delta = M_.params( 3 );
M_.params( 4 ) = 1;
lbar = M_.params( 4 );
M_.params( 5 ) = 2;
gamma = M_.params( 5 );
M_.params( 6 ) = 0;
mu = M_.params( 6 );
M_.params( 7 ) = 0.05;
sigma = M_.params( 7 );
M_.params( 8 ) = 0.9;
rho = M_.params( 8 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = exp(M_.params(6))*21.6333351205289^M_.params(1)*M_.params(4)^(1-M_.params(1));
oo_.steady_state( 2 ) = M_.params(4)*exp(M_.params(6))*(1-M_.params(1))*(21.6333351205289/M_.params(4))^M_.params(1)+21.6333351205289*(1+exp(M_.params(6))*M_.params(1)*(M_.params(4)/21.6333351205289)^(1-M_.params(1))-M_.params(3));
oo_.steady_state( 3 ) = 21.6333351205289;
oo_.steady_state( 5 ) = exp(M_.params(6))*M_.params(1)*(M_.params(4)/21.6333351205289)^(1-M_.params(1));
oo_.steady_state( 6 ) = exp(M_.params(6))*(1-M_.params(1))*(21.6333351205289/M_.params(4))^M_.params(1);
oo_.steady_state( 7 ) = M_.params(6);
oo_.endo_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
options_.solve_algo = 0;
steady;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(7))^2;
M_.sigma_e_is_diagonal = 1;
options_.irf = 200;
options_.order = 2;
var_list_=[];
info = stoch_simul(var_list_);
save TayApproxDyn.mat ;
save('TayApproxDyn_results.mat', 'oo_', 'M_', 'options_');
diary off

disp(['Total computing time : ' dynsec2hms(toc) ]);
