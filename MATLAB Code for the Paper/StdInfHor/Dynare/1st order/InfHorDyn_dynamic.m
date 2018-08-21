function [residual, g1, g2, g3] = InfHorDyn_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           columns: equations in order of declaration
%                                                           rows: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(6, 1);
T71 = params(5)^(1-params(1))*exp(y(6))*y(1)^params(1);
T75 = exp(y(6))*(1-params(1))*(y(1)/params(5))^params(1);
T79 = exp(y(6))*params(1)*(params(5)/y(1))^(1-params(1));
lhs =y(7);
rhs =y(3)*params(5)+(1+y(8)-params(3))*y(1)-y(5);
residual(1)= lhs-rhs;
lhs =y(7)^(-params(4));
rhs =params(2)*(1+y(10)-params(3))*y(9)^(-params(4));
residual(2)= lhs-rhs;
lhs =y(4);
rhs =T71;
residual(3)= lhs-rhs;
lhs =y(3);
rhs =T75;
residual(4)= lhs-rhs;
lhs =y(8);
rhs =T79;
residual(5)= lhs-rhs;
lhs =y(6);
rhs =params(7)*y(2)+(1-params(7))*params(6)+x(it_, 1);
residual(6)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(6, 11);

  %
  % Jacobian matrix
  %

  g1(1,3)=(-params(5));
  g1(1,1)=(-(1+y(8)-params(3)));
  g1(1,5)=1;
  g1(1,7)=1;
  g1(1,8)=(-y(1));
  g1(2,7)=getPowerDeriv(y(7),(-params(4)),1);
  g1(2,9)=(-(params(2)*(1+y(10)-params(3))*getPowerDeriv(y(9),(-params(4)),1)));
  g1(2,10)=(-(params(2)*y(9)^(-params(4))));
  g1(3,4)=1;
  g1(3,1)=(-(params(5)^(1-params(1))*exp(y(6))*getPowerDeriv(y(1),params(1),1)));
  g1(3,6)=(-T71);
  g1(4,3)=1;
  g1(4,1)=(-(exp(y(6))*(1-params(1))*1/params(5)*getPowerDeriv(y(1)/params(5),params(1),1)));
  g1(4,6)=(-T75);
  g1(5,1)=(-(exp(y(6))*params(1)*(-params(5))/(y(1)*y(1))*getPowerDeriv(params(5)/y(1),1-params(1),1)));
  g1(5,6)=(-T79);
  g1(5,8)=1;
  g1(6,2)=(-params(7));
  g1(6,6)=1;
  g1(6,11)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,121);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,1331);
end
end
