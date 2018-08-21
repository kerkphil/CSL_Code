%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'varyingMtax_d2';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('varyingMtax_d2.log');
M_.exo_names = 'eps';
M_.exo_names_tex = 'eps';
M_.exo_names_long = 'eps';
M_.endo_names = 'w';
M_.endo_names_tex = 'w';
M_.endo_names_long = 'w';
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'tau');
M_.endo_names_tex = char(M_.endo_names_tex, 'tau');
M_.endo_names_long = char(M_.endo_names_long, 'tau');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.param_names = 'calpha';
M_.param_names_tex = 'calpha';
M_.param_names_long = 'calpha';
M_.param_names = char(M_.param_names, 'cbeta');
M_.param_names_tex = char(M_.param_names_tex, 'cbeta');
M_.param_names_long = char(M_.param_names_long, 'cbeta');
M_.param_names = char(M_.param_names, 'cdelta');
M_.param_names_tex = char(M_.param_names_tex, 'cdelta');
M_.param_names_long = char(M_.param_names_long, 'cdelta');
M_.param_names = char(M_.param_names, 'cgamma');
M_.param_names_tex = char(M_.param_names_tex, 'cgamma');
M_.param_names_long = char(M_.param_names_long, 'cgamma');
M_.param_names = char(M_.param_names, 'cmu');
M_.param_names_tex = char(M_.param_names_tex, 'cmu');
M_.param_names_long = char(M_.param_names_long, 'cmu');
M_.param_names = char(M_.param_names, 'crho');
M_.param_names_tex = char(M_.param_names_tex, 'crho');
M_.param_names_long = char(M_.param_names_long, 'crho');
M_.param_names = char(M_.param_names, 'csigma');
M_.param_names_tex = char(M_.param_names_tex, 'csigma');
M_.param_names_long = char(M_.param_names_long, 'csigma');
M_.param_names = char(M_.param_names, 'ctau1');
M_.param_names_tex = char(M_.param_names_tex, 'ctau1');
M_.param_names_long = char(M_.param_names_long, 'ctau1');
M_.param_names = char(M_.param_names, 'ctau2');
M_.param_names_tex = char(M_.param_names_tex, 'ctau2');
M_.param_names_long = char(M_.param_names_long, 'ctau2');
M_.param_names = char(M_.param_names, 'ca');
M_.param_names_tex = char(M_.param_names_tex, 'ca');
M_.param_names_long = char(M_.param_names_long, 'ca');
M_.param_names = char(M_.param_names, 'cb');
M_.param_names_tex = char(M_.param_names_tex, 'cb');
M_.param_names_long = char(M_.param_names_long, 'cb');
M_.param_names = char(M_.param_names, 'cpi');
M_.param_names_tex = char(M_.param_names_tex, 'cpi');
M_.param_names_long = char(M_.param_names_long, 'cpi');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 7;
M_.param_nbr = 12;
M_.orig_endo_nbr = 7;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('varyingMtax_d2_static');
erase_compiled_function('varyingMtax_d2_dynamic');
M_.lead_lag_incidence = [
 0 3 10;
 1 4 0;
 2 5 0;
 0 6 11;
 0 7 12;
 0 8 13;
 0 9 0;]';
M_.nstatic = 1;
M_.nfwrd   = 4;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 4;
M_.nspred   = 2;
M_.ndynamic   = 6;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(12, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 31;
M_.NNZDerivatives(2) = 54;
M_.NNZDerivatives(3) = -1;
load '../ParamSet/varyingMtax_par.mat' alpha beta delta gamma mu rho sigma tau1 tau2 atau btau css rss wss tauss parpi ;
M_.params( 1 ) = alpha;
calpha = M_.params( 1 );
M_.params( 2 ) = beta;
cbeta = M_.params( 2 );
M_.params( 3 ) = delta;
cdelta = M_.params( 3 );
M_.params( 4 ) = gamma;
cgamma = M_.params( 4 );
M_.params( 5 ) = mu;
cmu = M_.params( 5 );
M_.params( 6 ) = rho;
crho = M_.params( 6 );
M_.params( 7 ) = sigma;
csigma = M_.params( 7 );
M_.params( 8 ) = tau1;
ctau1 = M_.params( 8 );
M_.params( 9 ) = tau2;
ctau2 = M_.params( 9 );
M_.params( 10 ) = atau;
ca = M_.params( 10 );
M_.params( 11 ) = btau;
cb = M_.params( 11 );
M_.params( 12 ) = parpi;
cpi = M_.params( 12 );
clear alpha beta delta gamma mu rho sigma T tau1 tau2 atau btau parpi ;
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 4 ) = 1.2961;
oo_.steady_state( 5 ) = 0.0971;
oo_.steady_state( 1 ) = 1.2961;
oo_.steady_state( 2 ) = 7.1838;
oo_.steady_state( 3 ) = M_.params(5);
oo_.steady_state( 6 ) = 0.1075;
oo_.steady_state( 7 ) = 0.1757;
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
steady;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(7))^2;
options_.irf = 50;
options_.order = 2;
var_list_=[];
info = stoch_simul(var_list_);
save('varyingMtax_d2_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('varyingMtax_d2_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('varyingMtax_d2_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('varyingMtax_d2_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('varyingMtax_d2_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
