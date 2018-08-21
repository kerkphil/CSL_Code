function [residual, g1, g2, g3] = StdInfHor_d2_dynamic(y, x, params, steady_state, it_)
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
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(5, 1);
lhs =y(6)^(-params(4));
rhs =params(2)*(1+y(9)-params(3))*y(8)^(-params(4));
residual(1)= lhs-rhs;
lhs =y(6);
rhs =y(3)+(1+y(7)-params(3))*y(1)-y(4);
residual(2)= lhs-rhs;
lhs =y(3);
rhs =(1-params(1))*exp(y(5))*y(1)^params(1);
residual(3)= lhs-rhs;
lhs =y(7);
rhs =params(1)*exp(y(5))*y(1)^(params(1)-1);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =params(6)*y(2)+(1-params(6))*params(5)+x(it_, 1);
residual(5)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(5, 10);

  %
  % Jacobian matrix
  %

  g1(1,6)=getPowerDeriv(y(6),(-params(4)),1);
  g1(1,8)=(-(params(2)*(1+y(9)-params(3))*getPowerDeriv(y(8),(-params(4)),1)));
  g1(1,9)=(-(params(2)*y(8)^(-params(4))));
  g1(2,3)=(-1);
  g1(2,1)=(-(1+y(7)-params(3)));
  g1(2,4)=1;
  g1(2,6)=1;
  g1(2,7)=(-y(1));
  g1(3,3)=1;
  g1(3,1)=(-((1-params(1))*exp(y(5))*getPowerDeriv(y(1),params(1),1)));
  g1(3,5)=(-((1-params(1))*exp(y(5))*y(1)^params(1)));
  g1(4,1)=(-(params(1)*exp(y(5))*getPowerDeriv(y(1),params(1)-1,1)));
  g1(4,5)=(-(params(1)*exp(y(5))*y(1)^(params(1)-1)));
  g1(4,7)=1;
  g1(5,2)=(-params(6));
  g1(5,5)=1;
  g1(5,10)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(14,3);
  v2(1,1)=1;
  v2(1,2)=56;
  v2(1,3)=getPowerDeriv(y(6),(-params(4)),2);
  v2(2,1)=1;
  v2(2,2)=78;
  v2(2,3)=(-(params(2)*(1+y(9)-params(3))*getPowerDeriv(y(8),(-params(4)),2)));
  v2(3,1)=1;
  v2(3,2)=88;
  v2(3,3)=(-(params(2)*getPowerDeriv(y(8),(-params(4)),1)));
  v2(4,1)=1;
  v2(4,2)=79;
  v2(4,3)=  v2(3,3);
  v2(5,1)=2;
  v2(5,2)=61;
  v2(5,3)=(-1);
  v2(6,1)=2;
  v2(6,2)=7;
  v2(6,3)=  v2(5,3);
  v2(7,1)=3;
  v2(7,2)=1;
  v2(7,3)=(-((1-params(1))*exp(y(5))*getPowerDeriv(y(1),params(1),2)));
  v2(8,1)=3;
  v2(8,2)=41;
  v2(8,3)=(-((1-params(1))*exp(y(5))*getPowerDeriv(y(1),params(1),1)));
  v2(9,1)=3;
  v2(9,2)=5;
  v2(9,3)=  v2(8,3);
  v2(10,1)=3;
  v2(10,2)=45;
  v2(10,3)=(-((1-params(1))*exp(y(5))*y(1)^params(1)));
  v2(11,1)=4;
  v2(11,2)=1;
  v2(11,3)=(-(params(1)*exp(y(5))*getPowerDeriv(y(1),params(1)-1,2)));
  v2(12,1)=4;
  v2(12,2)=41;
  v2(12,3)=(-(params(1)*exp(y(5))*getPowerDeriv(y(1),params(1)-1,1)));
  v2(13,1)=4;
  v2(13,2)=5;
  v2(13,3)=  v2(12,3);
  v2(14,1)=4;
  v2(14,2)=45;
  v2(14,3)=(-(params(1)*exp(y(5))*y(1)^(params(1)-1)));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),5,100);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],5,1000);
end
end
