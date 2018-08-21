function [residual, g1, g2] = StdInfHor_d1_static(y, x, params)
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
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 5, 1);

%
% Model equations
%

T11 = y(4)^(-params(4));
lhs =T11;
rhs =T11*params(2)*(1+y(5)-params(3));
residual(1)= lhs-rhs;
lhs =y(4);
rhs =y(1)+(1+y(5)-params(3))*y(2)-y(2);
residual(2)= lhs-rhs;
lhs =y(1);
rhs =(1-params(1))*exp(y(3))*y(2)^params(1);
residual(3)= lhs-rhs;
lhs =y(5);
rhs =params(1)*exp(y(3))*y(2)^(params(1)-1);
residual(4)= lhs-rhs;
lhs =y(3);
rhs =y(3)*params(6)+(1-params(6))*params(5)+x(1);
residual(5)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(5, 5);

  %
  % Jacobian matrix
  %

  g1(1,4)=getPowerDeriv(y(4),(-params(4)),1)-params(2)*(1+y(5)-params(3))*getPowerDeriv(y(4),(-params(4)),1);
  g1(1,5)=(-(T11*params(2)));
  g1(2,1)=(-1);
  g1(2,2)=(-(1+y(5)-params(3)-1));
  g1(2,4)=1;
  g1(2,5)=(-y(2));
  g1(3,1)=1;
  g1(3,2)=(-((1-params(1))*exp(y(3))*getPowerDeriv(y(2),params(1),1)));
  g1(3,3)=(-((1-params(1))*exp(y(3))*y(2)^params(1)));
  g1(4,2)=(-(params(1)*exp(y(3))*getPowerDeriv(y(2),params(1)-1,1)));
  g1(4,3)=(-(params(1)*exp(y(3))*y(2)^(params(1)-1)));
  g1(4,5)=1;
  g1(5,3)=1-params(6);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],5,25);
end
end
