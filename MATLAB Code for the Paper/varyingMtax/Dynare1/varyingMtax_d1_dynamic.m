function [residual, g1, g2, g3] = varyingMtax_d1_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(7, 1);
T15 = params(2)*y(11)^(-params(4));
T27 = (params(9)-params(8))/params(12);
T38 = 1+(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))^2;
T40 = T27*(y(12)-params(3))*params(10)/T38;
T75 = 1/params(12);
lhs =y(6)^(-params(4));
rhs =T15*(1+(1-y(13))*(y(12)-params(3))-(y(10)+(y(12)-params(3))*y(4))*T40);
residual(1)= lhs-rhs;
lhs =y(6);
rhs =y(3)+(1+y(7)-params(3))*y(1)-y(8)*(y(3)+y(1)*(y(7)-params(3)))-y(4)+y(9);
residual(2)= lhs-rhs;
lhs =y(3);
rhs =(1-params(1))*exp(y(5))*y(1)^params(1);
residual(3)= lhs-rhs;
lhs =y(7);
rhs =params(1)*exp(y(5))*y(1)^(params(1)-1);
residual(4)= lhs-rhs;
lhs =y(8);
rhs =params(8)+(params(9)-params(8))*(T75*atan(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))+0.5);
residual(5)= lhs-rhs;
lhs =y(9);
rhs =y(8)*(y(3)+y(1)*(y(7)-params(3)));
residual(6)= lhs-rhs;
lhs =y(5);
rhs =params(6)*y(2)+(1-params(6))*params(5)+x(it_, 1);
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 14);

  %
  % Jacobian matrix
  %

  g1(1,10)=(-(T15*(-(T40+(y(10)+(y(12)-params(3))*y(4))*T27*(-((y(12)-params(3))*params(10)*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))))/(T38*T38)))));
  g1(1,4)=(-(T15*(-((y(12)-params(3))*T40+(y(10)+(y(12)-params(3))*y(4))*T27*(-((y(12)-params(3))*params(10)*(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))))/(T38*T38)))));
  g1(1,6)=getPowerDeriv(y(6),(-params(4)),1);
  g1(1,11)=(-((1+(1-y(13))*(y(12)-params(3))-(y(10)+(y(12)-params(3))*y(4))*T40)*params(2)*getPowerDeriv(y(11),(-params(4)),1)));
  g1(1,12)=(-(T15*(1-y(13)-(y(4)*T40+(y(10)+(y(12)-params(3))*y(4))*T27*(params(10)*T38-(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))*params(10)*y(4))/(T38*T38)))));
  g1(1,13)=(-(T15*(-(y(12)-params(3)))));
  g1(2,3)=(-(1-y(8)));
  g1(2,1)=(-(1+y(7)-params(3)-y(8)*(y(7)-params(3))));
  g1(2,4)=1;
  g1(2,6)=1;
  g1(2,7)=(-(y(1)-y(1)*y(8)));
  g1(2,8)=y(3)+y(1)*(y(7)-params(3));
  g1(2,9)=(-1);
  g1(3,3)=1;
  g1(3,1)=(-((1-params(1))*exp(y(5))*getPowerDeriv(y(1),params(1),1)));
  g1(3,5)=(-((1-params(1))*exp(y(5))*y(1)^params(1)));
  g1(4,1)=(-(params(1)*exp(y(5))*getPowerDeriv(y(1),params(1)-1,1)));
  g1(4,5)=(-(params(1)*exp(y(5))*y(1)^(params(1)-1)));
  g1(4,7)=1;
  g1(5,3)=(-((params(9)-params(8))*T75*params(10)/(1+(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3)))))));
  g1(5,1)=(-((params(9)-params(8))*T75*params(10)*(y(7)-params(3))/(1+(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3)))))));
  g1(5,7)=(-((params(9)-params(8))*T75*params(10)*y(1)/(1+(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3)))))));
  g1(5,8)=1;
  g1(6,3)=(-y(8));
  g1(6,1)=(-(y(8)*(y(7)-params(3))));
  g1(6,7)=(-(y(1)*y(8)));
  g1(6,8)=(-(y(3)+y(1)*(y(7)-params(3))));
  g1(6,9)=1;
  g1(7,2)=(-params(6));
  g1(7,5)=1;
  g1(7,14)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,196);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,2744);
end
end
