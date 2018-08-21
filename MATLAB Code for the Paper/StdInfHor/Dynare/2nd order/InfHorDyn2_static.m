function [residual, g1, g2] = InfHorDyn2_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                    columns: equations in order of declaration
%                                                    rows: variables in declaration order
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: equations in order of declaration
%                                                       rows: variables in declaration order
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 6, 1);

%
% Model equations
%

T23 = y(5)^(-params(4));
T36 = params(5)^(1-params(1))*exp(y(4))*y(3)^params(1);
T41 = (1-params(1))*exp(y(4))*(y(3)/params(5))^params(1);
T46 = params(1)*exp(y(4))*(params(5)/y(3))^(1-params(1));
lhs =y(5);
rhs =y(1)*params(5)+(1+y(6)-params(3))*y(3)-y(3);
residual(1)= lhs-rhs;
lhs =T23;
rhs =params(2)*(1+y(6)-params(3))*T23;
residual(2)= lhs-rhs;
lhs =y(2);
rhs =T36;
residual(3)= lhs-rhs;
lhs =y(1);
rhs =T41;
residual(4)= lhs-rhs;
lhs =y(6);
rhs =T46;
residual(5)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(6)+x(1);
residual(6)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(6, 6);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-params(5));
  g1(1,3)=(-(1+y(6)-params(3)-1));
  g1(1,5)=1;
  g1(1,6)=(-y(3));
  g1(2,5)=getPowerDeriv(y(5),(-params(4)),1)-params(2)*(1+y(6)-params(3))*getPowerDeriv(y(5),(-params(4)),1);
  g1(2,6)=(-(T23*params(2)));
  g1(3,2)=1;
  g1(3,3)=(-(params(5)^(1-params(1))*exp(y(4))*getPowerDeriv(y(3),params(1),1)));
  g1(3,4)=(-T36);
  g1(4,1)=1;
  g1(4,3)=(-((1-params(1))*exp(y(4))*1/params(5)*getPowerDeriv(y(3)/params(5),params(1),1)));
  g1(4,4)=(-T41);
  g1(5,3)=(-(params(1)*exp(y(4))*(-params(5))/(y(3)*y(3))*getPowerDeriv(params(5)/y(3),1-params(1),1)));
  g1(5,4)=(-T46);
  g1(5,6)=1;
  g1(6,4)=1-params(6);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,36);
end
end
