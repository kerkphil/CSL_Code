function ee = CSL1ss(Kbar)
% returns the value of the Euler equation for the input value of Kbar
% used in conjunction with fsolve in CSL.m to find the value of Kbar

in = [Kbar; Kbar; Kbar; 0; 0];

ee = CSLStdInfHordyn(in);