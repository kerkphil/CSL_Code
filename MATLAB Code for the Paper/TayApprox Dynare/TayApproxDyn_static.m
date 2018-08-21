function [residual, g1, g2] = TayApproxDyn_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 7, 1);

%
% Model equations
%

T11 = y(2)^(-params(5));
T36 = exp(y(7))*y(3)^params(1)*params(4)^(1-params(1));
T44 = exp(y(7))*params(1)*(params(4)/y(3))^(1-params(1));
T49 = exp(y(7))*(1-params(1))*(y(3)/params(4))^params(1);
lhs =T11;
rhs =T11*params(2)*(1+y(5)-params(3));
residual(1)= lhs-rhs;
lhs =y(2)+y(3);
rhs =y(6)*params(4)+(1+y(5)-params(3))*y(3);
residual(2)= lhs-rhs;
lhs =y(1);
rhs =T36;
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(1)-y(2);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =T44;
residual(5)= lhs-rhs;
lhs =y(6);
rhs =T49;
residual(6)= lhs-rhs;
lhs =y(7);
rhs =y(7)*params(8)+(1-params(8))*params(6)+x(1);
residual(7)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(7, 7);

%
% Jacobian matrix
%

  g1(1,2)=getPowerDeriv(y(2),(-params(5)),1)-params(2)*(1+y(5)-params(3))*getPowerDeriv(y(2),(-params(5)),1);
  g1(1,5)=(-(T11*params(2)));
  g1(2,2)=1;
  g1(2,3)=1-(1+y(5)-params(3));
  g1(2,5)=(-y(3));
  g1(2,6)=(-params(4));
  g1(3,1)=1;
  g1(3,3)=(-(params(4)^(1-params(1))*exp(y(7))*getPowerDeriv(y(3),params(1),1)));
  g1(3,7)=(-T36);
  g1(4,1)=(-1);
  g1(4,2)=1;
  g1(4,4)=1;
  g1(5,3)=(-(exp(y(7))*params(1)*(-params(4))/(y(3)*y(3))*getPowerDeriv(params(4)/y(3),1-params(1),1)));
  g1(5,5)=1;
  g1(5,7)=(-T44);
  g1(6,3)=(-(exp(y(7))*(1-params(1))*1/params(4)*getPowerDeriv(y(3)/params(4),params(1),1)));
  g1(6,6)=1;
  g1(6,7)=(-T49);
  g1(7,7)=1-params(8);
  if ~isreal(g1)
    g1 = real(g1)+imag(g1).^2;
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],7,49);
end
end
