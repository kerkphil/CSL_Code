% Solving the stochastic neoclassical growth model by hand:

disp('EXAMPLE0: The stochastic neoclassical growth model, solved by hand.');

% Setting parameters:

CHOICE = 0;
if CHOICE,
   n = 1;
   m = 1;
else
   n = 4;  
   m = 5;
end;

o = ones(n,m);
o_m = ones(1,m);
o_n = ones(n,1);

Z_bar     = 1*o;       % Normalization
rho       = .36*o;     % Capital share
if CHOICE,
   delta     = .025*o; % Depreciation rate for capital
   eta       = 1.0*o;  % constant of relative risk aversion 
                       % = 1/(coeff. of intertemporal substitution)
else
   eta = o_n*[ .01, .5, 1, 2, 1000 ];
   delta = [ 0.0, 0.025, 0.1, 1.0 ]'*o_m;
end;   
R_bar     = 1.01*o;    % One percent real interest per quarter
psi       = .95*o;     % autocorrelation of technology shock
sigma_eps = .712*o;    % Standard deviation of technology shock.  Units: Percent.

% Calculating the steady state:

betta   = 1.0./R_bar;  % Discount factor beta
K_bar   = ((rho.*Z_bar)./(R_bar - 1 + delta)).^(1.0./(1 - rho));
Y_bar   = Z_bar.*K_bar.^rho;
C_bar   = Y_bar - delta.*K_bar;

gamma = (1-betta.*(1-delta)).*(1-rho).*(1-betta+betta.*delta.*(1-rho))./(eta.*rho.*betta) ...
         + 1 + 1.0./betta;
nu_kk = gamma/2 - ( (gamma/2).^2 - 1.0./betta).^(.5);
nu_ck = (1.0./betta - nu_kk).*K_bar./C_bar;
nu_rk = - (1 - betta.*(1-delta)).*(1 - rho);

nu_rz = 1 - betta.*(1 - delta);
nu_kz = (nu_rz.*psi + eta.*(1-psi).*Y_bar./C_bar)./(-nu_rk+eta.*nu_ck + eta.*(1-psi).*K_bar./C_bar);
nu_cz = Y_bar./C_bar - K_bar./C_bar .* nu_kz;


if CHOICE,
   disp(sprintf('nu_kk=%6.3f, nu_kz=%6.3f',nu_kk,nu_kz));
   disp(sprintf('nu_ck=%6.3f, nu_cz=%6.3f',nu_ck,nu_cz));
   disp(sprintf('nu_rk=%6.3f, nu_rz=%6.3f',nu_rk,nu_rz));
else
   nu_kk
   nu_kz
end;


