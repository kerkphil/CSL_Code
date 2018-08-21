function [residual, g1, g2, g3] = varyingMtax_d2_dynamic(y, x, params, steady_state, it_)
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
T99 = 1+(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))));
T111 = T27*(-((y(12)-params(3))*params(10)*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))))/(T38*T38);
T136 = T27*(-((y(12)-params(3))*params(10)*(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))))/(T38*T38);
T148 = params(2)*getPowerDeriv(y(11),(-params(4)),1);
T164 = params(10)*T38-(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))*params(10)*y(4);
T166 = T27*T164/(T38*T38);
T182 = T38*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))+T38*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11));
T216 = T38*(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))+T38*(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11));
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

  g1(1,10)=(-(T15*(-(T40+(y(10)+(y(12)-params(3))*y(4))*T111))));
  g1(1,4)=(-(T15*(-((y(12)-params(3))*T40+(y(10)+(y(12)-params(3))*y(4))*T136))));
  g1(1,6)=getPowerDeriv(y(6),(-params(4)),1);
  g1(1,11)=(-((1+(1-y(13))*(y(12)-params(3))-(y(10)+(y(12)-params(3))*y(4))*T40)*T148));
  g1(1,12)=(-(T15*(1-y(13)-(y(4)*T40+(y(10)+(y(12)-params(3))*y(4))*T166))));
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
  g1(5,3)=(-((params(9)-params(8))*T75*params(10)/T99));
  g1(5,1)=(-((params(9)-params(8))*T75*params(10)*(y(7)-params(3))/T99));
  g1(5,7)=(-((params(9)-params(8))*T75*params(10)*y(1)/T99));
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

  v2 = zeros(54,3);
  v2(1,1)=1;
  v2(1,2)=136;
  v2(1,3)=(-(T15*(-(T111+T111+(y(10)+(y(12)-params(3))*y(4))*T27*(T38*T38*(-((y(12)-params(3))*params(10)*params(10)*2*params(10)))-(-((y(12)-params(3))*params(10)*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))))*T182)/(T38*T38*T38*T38)))));
  v2(2,1)=1;
  v2(2,2)=52;
  v2(2,3)=(-(T15*(-((y(12)-params(3))*T111+T136+(y(10)+(y(12)-params(3))*y(4))*T27*(T38*T38*(-((y(12)-params(3))*params(10)*(y(12)-params(3))*params(10)*2*params(10)))-(-((y(12)-params(3))*params(10)*(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))))*T182)/(T38*T38*T38*T38)))));
  v2(3,1)=1;
  v2(3,2)=130;
  v2(3,3)=  v2(2,3);
  v2(4,1)=1;
  v2(4,2)=46;
  v2(4,3)=(-(T15*(-((y(12)-params(3))*T136+(y(12)-params(3))*T136+(y(10)+(y(12)-params(3))*y(4))*T27*(T38*T38*(-((y(12)-params(3))*params(10)*(y(12)-params(3))*params(10)*2*(y(12)-params(3))*params(10)))-(-((y(12)-params(3))*params(10)*(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))))*T216)/(T38*T38*T38*T38)))));
  v2(5,1)=1;
  v2(5,2)=76;
  v2(5,3)=getPowerDeriv(y(6),(-params(4)),2);
  v2(6,1)=1;
  v2(6,2)=150;
  v2(6,3)=(-((-(T40+(y(10)+(y(12)-params(3))*y(4))*T111))*T148));
  v2(7,1)=1;
  v2(7,2)=137;
  v2(7,3)=  v2(6,3);
  v2(8,1)=1;
  v2(8,2)=144;
  v2(8,3)=(-((-((y(12)-params(3))*T40+(y(10)+(y(12)-params(3))*y(4))*T136))*T148));
  v2(9,1)=1;
  v2(9,2)=53;
  v2(9,3)=  v2(8,3);
  v2(10,1)=1;
  v2(10,2)=151;
  v2(10,3)=(-((1+(1-y(13))*(y(12)-params(3))-(y(10)+(y(12)-params(3))*y(4))*T40)*params(2)*getPowerDeriv(y(11),(-params(4)),2)));
  v2(11,1)=1;
  v2(11,2)=164;
  v2(11,3)=(-(T15*(-(y(4)*T111+T166+(y(10)+(y(12)-params(3))*y(4))*T27*(T38*T38*(params(10)*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))-(y(12)-params(3))*params(10)*params(10)*y(4)*2*params(10))-T164*T182)/(T38*T38*T38*T38)))));
  v2(12,1)=1;
  v2(12,2)=138;
  v2(12,3)=  v2(11,3);
  v2(13,1)=1;
  v2(13,2)=158;
  v2(13,3)=(-(T15*(-(T40+y(4)*T136+(y(12)-params(3))*T166+(y(10)+(y(12)-params(3))*y(4))*T27*(T38*T38*(params(10)*(y(12)-params(3))*params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))-(y(12)-params(3))*params(10)*(params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))+params(10)*y(4)*2*(y(12)-params(3))*params(10)))-T164*T216)/(T38*T38*T38*T38)))));
  v2(14,1)=1;
  v2(14,2)=54;
  v2(14,3)=  v2(13,3);
  v2(15,1)=1;
  v2(15,2)=165;
  v2(15,3)=(-(T148*(1-y(13)-(y(4)*T40+(y(10)+(y(12)-params(3))*y(4))*T166))));
  v2(16,1)=1;
  v2(16,2)=152;
  v2(16,3)=  v2(15,3);
  v2(17,1)=1;
  v2(17,2)=166;
  v2(17,3)=(-(T15*(-(y(4)*T166+y(4)*T166+(y(10)+(y(12)-params(3))*y(4))*T27*(T38*T38*(params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))*params(10)*y(4)-(params(10)*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))*params(10)*y(4)+(y(12)-params(3))*params(10)*params(10)*y(4)*2*params(10)*y(4)))-T164*(T38*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))*params(10)*y(4)+T38*2*(params(10)*(y(10)+(y(12)-params(3))*y(4))+params(11))*params(10)*y(4)))/(T38*T38*T38*T38)))));
  v2(18,1)=1;
  v2(18,2)=179;
  v2(18,3)=(-(T148*(-(y(12)-params(3)))));
  v2(19,1)=1;
  v2(19,2)=153;
  v2(19,3)=  v2(18,3);
  v2(20,1)=1;
  v2(20,2)=180;
  v2(20,3)=T15;
  v2(21,1)=1;
  v2(21,2)=167;
  v2(21,3)=  v2(20,3);
  v2(22,1)=2;
  v2(22,2)=85;
  v2(22,3)=(-(1-y(8)));
  v2(23,1)=2;
  v2(23,2)=7;
  v2(23,3)=  v2(22,3);
  v2(24,1)=2;
  v2(24,2)=101;
  v2(24,3)=1;
  v2(25,1)=2;
  v2(25,2)=36;
  v2(25,3)=  v2(24,3);
  v2(26,1)=2;
  v2(26,2)=99;
  v2(26,3)=y(7)-params(3);
  v2(27,1)=2;
  v2(27,2)=8;
  v2(27,3)=  v2(26,3);
  v2(28,1)=2;
  v2(28,2)=105;
  v2(28,3)=y(1);
  v2(29,1)=2;
  v2(29,2)=92;
  v2(29,3)=  v2(28,3);
  v2(30,1)=3;
  v2(30,2)=1;
  v2(30,3)=(-((1-params(1))*exp(y(5))*getPowerDeriv(y(1),params(1),2)));
  v2(31,1)=3;
  v2(31,2)=57;
  v2(31,3)=(-((1-params(1))*exp(y(5))*getPowerDeriv(y(1),params(1),1)));
  v2(32,1)=3;
  v2(32,2)=5;
  v2(32,3)=  v2(31,3);
  v2(33,1)=3;
  v2(33,2)=61;
  v2(33,3)=(-((1-params(1))*exp(y(5))*y(1)^params(1)));
  v2(34,1)=4;
  v2(34,2)=1;
  v2(34,3)=(-(params(1)*exp(y(5))*getPowerDeriv(y(1),params(1)-1,2)));
  v2(35,1)=4;
  v2(35,2)=57;
  v2(35,3)=(-(params(1)*exp(y(5))*getPowerDeriv(y(1),params(1)-1,1)));
  v2(36,1)=4;
  v2(36,2)=5;
  v2(36,3)=  v2(35,3);
  v2(37,1)=4;
  v2(37,2)=61;
  v2(37,3)=(-(params(1)*exp(y(5))*y(1)^(params(1)-1)));
  v2(38,1)=5;
  v2(38,2)=31;
  v2(38,3)=(-((params(9)-params(8))*T75*(-(params(10)*(params(10)*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))+params(10)*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3)))))))/(T99*T99)));
  v2(39,1)=5;
  v2(39,2)=3;
  v2(39,3)=(-((params(9)-params(8))*T75*(-(params(10)*(y(7)-params(3))*(params(10)*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))+params(10)*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3)))))))/(T99*T99)));
  v2(40,1)=5;
  v2(40,2)=29;
  v2(40,3)=  v2(39,3);
  v2(41,1)=5;
  v2(41,2)=1;
  v2(41,3)=(-((params(9)-params(8))*T75*(-(params(10)*(y(7)-params(3))*((params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*params(10)*(y(7)-params(3))+(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*params(10)*(y(7)-params(3)))))/(T99*T99)));
  v2(42,1)=5;
  v2(42,2)=87;
  v2(42,3)=(-((params(9)-params(8))*T75*(-(params(10)*y(1)*(params(10)*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))+params(10)*(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3)))))))/(T99*T99)));
  v2(43,1)=5;
  v2(43,2)=35;
  v2(43,3)=  v2(42,3);
  v2(44,1)=5;
  v2(44,2)=85;
  v2(44,3)=(-((params(9)-params(8))*T75*(params(10)*T99-params(10)*y(1)*((params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*params(10)*(y(7)-params(3))+(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*params(10)*(y(7)-params(3))))/(T99*T99)));
  v2(45,1)=5;
  v2(45,2)=7;
  v2(45,3)=  v2(44,3);
  v2(46,1)=5;
  v2(46,2)=91;
  v2(46,3)=(-((params(9)-params(8))*T75*(-(params(10)*y(1)*((params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*params(10)*y(1)+(params(11)+params(10)*(y(3)+y(1)*(y(7)-params(3))))*params(10)*y(1))))/(T99*T99)));
  v2(47,1)=6;
  v2(47,2)=85;
  v2(47,3)=(-y(8));
  v2(48,1)=6;
  v2(48,2)=7;
  v2(48,3)=  v2(47,3);
  v2(49,1)=6;
  v2(49,2)=101;
  v2(49,3)=(-1);
  v2(50,1)=6;
  v2(50,2)=36;
  v2(50,3)=  v2(49,3);
  v2(51,1)=6;
  v2(51,2)=99;
  v2(51,3)=(-(y(7)-params(3)));
  v2(52,1)=6;
  v2(52,2)=8;
  v2(52,3)=  v2(51,3);
  v2(53,1)=6;
  v2(53,2)=105;
  v2(53,3)=(-y(1));
  v2(54,1)=6;
  v2(54,2)=92;
  v2(54,3)=  v2(53,3);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),7,196);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,2744);
end
end
