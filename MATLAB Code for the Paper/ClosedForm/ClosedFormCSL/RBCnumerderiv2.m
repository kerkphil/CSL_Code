function [AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, TT] = ...
    RBCnumerderiv2(funcname,theta0,params)

%This function takes the following inputs:
%  funcname - is the name of the function which generates a column vector
%  from ny+nx dynamic equations 
%    with the ny eqyations to be linearized into AX(t)+BX(t-1)+CY(t)+DZ(t)=0 in the first ny rows
%    the function must be written so as to evaluate to zero for all rows in the steady state
%  xybar - is a matrix of nx+ny steady state values 
%    with the values of xbar in the first nx rows and the values of ybar in the last ny rows
%  nx - is the number of elements in X
%  ny - is the number of elements in Y
%  nz - is the number of elements in Z
%The function generates matrices AA thru MM from the lineaization of the equations in the function "funcname"
%Note this subroutine does NOT find log-linearized derivatives
nx   = params(6) ;
ny   = params(7) ;
nz   = params(8) ;
incr = params(9) ;    %increment to use in calculating numerical derivatives

T0 = funcname(theta0,params);  %this should be zero or very close to it when evaluating at the SS

% create matrix of deviations for input, dev 
% deviate columns by adding incr
leng = 3*nx+2*(ny+nz);
dev=repmat(theta0',leng,1);
for i=1:leng
  dev(i,i)=dev(i,i)+incr;
end

% calculate a matrix of deviations for dynamic equations, big
% rows correspond to equations
% colums correspond to variables from "theta0" vector being changed
% note output of RBC2dyn is a column vector
big = zeros(nx+ny,leng);
% note the functions evaluate to 0 in the SS, not 1 as they would with log-linearization.
for i = 2:leng
    big(:,i) = funcname(dev(i,:)',params)-T0;
end
big = big/incr;

% pull out appropriate parts of "big" into Uhlig's matrices
AA = big(1:ny,nx+1:2*nx);
BB = big(1:ny,2*nx+1:3*nx);
CC = big(1:ny,3*nx+ny+1:3*nx+2*ny);
DD = big(1:ny,3*nx+2*ny+nz+1:leng);
FF = big(ny+1:ny+nx,1:nx);
GG = big(ny+1:ny+nx,nx+1:2*nx);
HH = big(ny+1:ny+nx,2*nx+1:3*nx);
JJ = big(ny+1:ny+nx,3*nx+1:3*nx+ny);
KK = big(ny+1:ny+nx,3*nx+ny+1:3*nx+2*ny);
LL = big(ny+1:ny+nx,3*nx+2*ny+1:3*nx+2*ny+nz);
MM = big(ny+1:ny+nx,3*nx+2*ny+nz+1:leng);
TT = T0;

end