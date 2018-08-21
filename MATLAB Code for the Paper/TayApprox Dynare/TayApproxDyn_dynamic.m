function [residual, g1, g2, g3] = TayApproxDyn_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(7, 1);
T41 = params(4)^(1-params(1));
T42 = exp(y(9))*y(1)^params(1)*T41;
T50 = exp(y(9))*params(1)*(params(4)/y(1))^(1-params(1));
T55 = exp(y(9))*(1-params(1))*(y(1)/params(4))^params(1);
T79 = getPowerDeriv(params(4)/y(1),1-params(1),1);
T87 = (-(exp(y(9))*(1-params(1))*1/params(4)*getPowerDeriv(y(1)/params(4),params(1),1)));
lhs =y(4)^(-params(5));
rhs =params(2)*(1+y(11)-params(3))*y(10)^(-params(5));
residual(1)= lhs-rhs;
lhs =y(4)+y(5);
rhs =y(8)*params(4)+(1+y(7)-params(3))*y(1);
residual(2)= lhs-rhs;
lhs =y(3);
rhs =T42;
residual(3)= lhs-rhs;
lhs =y(6);
rhs =y(3)-y(4);
residual(4)= lhs-rhs;
lhs =y(7);
rhs =T50;
residual(5)= lhs-rhs;
lhs =y(8);
rhs =T55;
residual(6)= lhs-rhs;
lhs =y(9);
rhs =params(8)*y(2)+(1-params(8))*params(6)+x(it_, 1);
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 12);

%
% Jacobian matrix
%

g1(1,4)=getPowerDeriv(y(4),(-params(5)),1);
g1(1,10)=(-(params(2)*(1+y(11)-params(3))*getPowerDeriv(y(10),(-params(5)),1)));
g1(1,11)=(-(params(2)*y(10)^(-params(5))));
g1(2,4)=1;
g1(2,1)=(-(1+y(7)-params(3)));
g1(2,5)=1;
g1(2,7)=(-y(1));
g1(2,8)=(-params(4));
g1(3,3)=1;
g1(3,1)=(-(T41*exp(y(9))*getPowerDeriv(y(1),params(1),1)));
g1(3,9)=(-T42);
g1(4,3)=(-1);
g1(4,4)=1;
g1(4,6)=1;
g1(5,1)=(-(exp(y(9))*params(1)*(-params(4))/(y(1)*y(1))*T79));
g1(5,7)=1;
g1(5,9)=(-T50);
g1(6,1)=T87;
g1(6,8)=1;
g1(6,9)=(-T55);
g1(7,2)=(-params(8));
g1(7,9)=1;
g1(7,12)=(-1);
end
if nargout >= 3,
%
% Hessian matrix
%

  v2 = zeros(18,3);
v2(1,1)=1;
v2(1,2)=40;
v2(1,3)=getPowerDeriv(y(4),(-params(5)),2);
v2(2,1)=1;
v2(2,2)=118;
v2(2,3)=(-(params(2)*(1+y(11)-params(3))*getPowerDeriv(y(10),(-params(5)),2)));
v2(3,1)=1;
v2(3,2)=130;
v2(3,3)=(-(params(2)*getPowerDeriv(y(10),(-params(5)),1)));
v2(4,1)=1;
v2(4,2)=119;
v2(4,3)=v2(3,3);
v2(5,1)=2;
v2(5,2)=73;
v2(5,3)=(-1);
v2(6,1)=2;
v2(6,2)=7;
v2(6,3)=v2(5,3);
v2(7,1)=3;
v2(7,2)=1;
v2(7,3)=(-(T41*exp(y(9))*getPowerDeriv(y(1),params(1),2)));
v2(8,1)=3;
v2(8,2)=97;
v2(8,3)=(-(T41*exp(y(9))*getPowerDeriv(y(1),params(1),1)));
v2(9,1)=3;
v2(9,2)=9;
v2(9,3)=v2(8,3);
v2(10,1)=3;
v2(10,2)=105;
v2(10,3)=(-T42);
v2(11,1)=5;
v2(11,2)=1;
v2(11,3)=(-(exp(y(9))*params(1)*(T79*(-((-params(4))*(y(1)+y(1))))/(y(1)*y(1)*y(1)*y(1))+(-params(4))/(y(1)*y(1))*(-params(4))/(y(1)*y(1))*getPowerDeriv(params(4)/y(1),1-params(1),2))));
v2(12,1)=5;
v2(12,2)=97;
v2(12,3)=(-(exp(y(9))*params(1)*(-params(4))/(y(1)*y(1))*T79));
v2(13,1)=5;
v2(13,2)=9;
v2(13,3)=v2(12,3);
v2(14,1)=5;
v2(14,2)=105;
v2(14,3)=(-T50);
v2(15,1)=6;
v2(15,2)=1;
v2(15,3)=(-(exp(y(9))*(1-params(1))*1/params(4)*1/params(4)*getPowerDeriv(y(1)/params(4),params(1),2)));
v2(16,1)=6;
v2(16,2)=97;
v2(16,3)=T87;
v2(17,1)=6;
v2(17,2)=9;
v2(17,3)=v2(16,3);
v2(18,1)=6;
v2(18,2)=105;
v2(18,3)=(-T55);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),7,144);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],7,1728);
end
end
