function ee = Mod1ss(Kbar,params)
% returns the value of the Euler equation for the input value of Kbar
% used in conjunction with fsolve in CSL.m to find the value of Kbar

in = [Kbar; Kbar; Kbar; 0; 0];

ee = Mod1dyn(in,params);