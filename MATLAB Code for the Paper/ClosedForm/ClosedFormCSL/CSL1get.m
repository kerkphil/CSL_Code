function [X,PP,QQ,UU] = CSL1get(X0,Z0,params)
%-------------------------------------------------------------------------%
% Calculates the value of X(t) given input values for X(t-1) & z(t) using
% a log-linear or simple linear approximation of the transition function
% about the current state X(t-1), z(t).
%-------------------------------------------------------------------------%

% find approximate linear policy function by linearizing about a given
% point
in = [X0; X0; X0; Z0; Z0];
[PP,QQ,UU] = CSL1findPQU(in,params);
X = X0.*exp(UU');  %log-linearized
% X = X0 + UU';  %linearized
