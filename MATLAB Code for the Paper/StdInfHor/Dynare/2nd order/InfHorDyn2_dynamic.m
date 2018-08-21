function [residual, g1, g2, g3] = InfHorDyn2_dynamic(y, x, params, steady_state, it_)
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
T41 = params(5)^(1-params(1));
T67 = T41*exp(y(6))*y(1)^params(1);
T71 = exp(y(6))*(1-params(1))*(y(1)/params(5))^params(1);
T75 = exp(y(6))*params(1)*(params(5)/y(1))^(1-params(1));
T87 = (-(exp(y(6))*(1-params(1))*1/params(5)*getPowerDeriv(y(1)/params(5),params(1),1)));
T90 = getPowerDeriv(params(5)/y(1),1-params(1),1);
lhs =y(7);
rhs =y(3)*params(5)+(1+y(8)-params(3))*y(1)-y(5);
residual(1)= lhs-rhs;
lhs =y(7)^(-params(4));
rhs =params(2)*y(9)^(-params(4))*(1+y(10)-params(3));
residual(2)= lhs-rhs;
lhs =y(4);
rhs =T67;
residual(3)= lhs-rhs;
lhs =y(3);
rhs =T71;
residual(4)= lhs-rhs;
lhs =y(8);
rhs =T75;
residual(5)= lhs-rhs;
lhs =y(6);
rhs =params(6)*y(2)+x(it_, 1);
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
  g1(3,1)=(-(T41*exp(y(6))*getPowerDeriv(y(1),params(1),1)));
  g1(3,6)=(-T67);
  g1(4,3)=1;
  g1(4,1)=T87;
  g1(4,6)=(-T71);
  g1(5,1)=(-(exp(y(6))*params(1)*(-params(5))/(y(1)*y(1))*T90));
  g1(5,6)=(-T75);
  g1(5,8)=1;
  g1(6,2)=(-params(6));
  g1(6,6)=1;
  g1(6,11)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(18,3);
  v2(1,1)=1;
  v2(1,2)=78;
  v2(1,3)=(-1);
  v2(2,1)=1;
  v2(2,2)=8;
  v2(2,3)=  v2(1,3);
  v2(3,1)=2;
  v2(3,2)=73;
  v2(3,3)=getPowerDeriv(y(7),(-params(4)),2);
  v2(4,1)=2;
  v2(4,2)=97;
  v2(4,3)=(-(params(2)*(1+y(10)-params(3))*getPowerDeriv(y(9),(-params(4)),2)));
  v2(5,1)=2;
  v2(5,2)=108;
  v2(5,3)=(-(params(2)*getPowerDeriv(y(9),(-params(4)),1)));
  v2(6,1)=2;
  v2(6,2)=98;
  v2(6,3)=  v2(5,3);
  v2(7,1)=3;
  v2(7,2)=1;
  v2(7,3)=(-(T41*exp(y(6))*getPowerDeriv(y(1),params(1),2)));
  v2(8,1)=3;
  v2(8,2)=56;
  v2(8,3)=(-(T41*exp(y(6))*getPowerDeriv(y(1),params(1),1)));
  v2(9,1)=3;
  v2(9,2)=6;
  v2(9,3)=  v2(8,3);
  v2(10,1)=3;
  v2(10,2)=61;
  v2(10,3)=(-T67);
  v2(11,1)=4;
  v2(11,2)=1;
  v2(11,3)=(-(exp(y(6))*(1-params(1))*1/params(5)*1/params(5)*getPowerDeriv(y(1)/params(5),params(1),2)));
  v2(12,1)=4;
  v2(12,2)=56;
  v2(12,3)=T87;
  v2(13,1)=4;
  v2(13,2)=6;
  v2(13,3)=  v2(12,3);
  v2(14,1)=4;
  v2(14,2)=61;
  v2(14,3)=(-T71);
  v2(15,1)=5;
  v2(15,2)=1;
  v2(15,3)=(-(exp(y(6))*params(1)*(T90*(-((-params(5))*(y(1)+y(1))))/(y(1)*y(1)*y(1)*y(1))+(-params(5))/(y(1)*y(1))*(-params(5))/(y(1)*y(1))*getPowerDeriv(params(5)/y(1),1-params(1),2))));
  v2(16,1)=5;
  v2(16,2)=56;
  v2(16,3)=(-(exp(y(6))*params(1)*(-params(5))/(y(1)*y(1))*T90));
  v2(17,1)=5;
  v2(17,2)=6;
  v2(17,3)=  v2(16,3);
  v2(18,1)=5;
  v2(18,2)=61;
  v2(18,3)=(-T75);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),6,121);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,1331);
end
end
