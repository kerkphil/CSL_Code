%-------------------------------------------------------------------------%
% Plot IRFs for one standard deviation shock for 4 different solutions
% methods of standard infinite horizon model
%-------------------------------------------------------------------------%
% This m-file calls the following additional m-file(s) and data file(s)
%   * 
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
cd('/Users/cc7768/Dropbox/Current State Linearization/MatLab/Std Inf Hor/Plotting') ;

%-------------------------------------------------------------------------%
% Pull in 4 series
%-------------------------------------------------------------------------%
cd('C:\Users\cc7768\Documents\GitHub\csl\MatLab\StdInfHor\Dynare\1st order') ;
load 'kmatirf1ord.mat' K_eps kss ;
K_eps_dyn1 = [0; 0; K_eps(1:end-1)] + kss * ones(length(K_eps)+1,1) ;
clear K_eps ;

cd('C:\Users\cc7768\Documents\GitHub\csl\MatLab\StdInfHor\Dynare\2nd order') ;
load 'kmatirf2ord.mat' K_eps kss2;
K_eps_dyn2 = [0; 0; K_eps(1:end-1)] + kss * ones(length(K_eps)+1,1) ;
clear K_eps ;

cd('C:\Users\cc7768\Documents\GitHub\csl\MatLab\StdInfHor\CSL') ;
load 'kmatirfcsl.mat' Kcs ;
K_eps_csl = Kcs ;
clear Kcs ;

cd('C:\Users\cc7768\Documents\GitHub\csl\MatLab\StdInfHor\Plotting') ;

T = length(K_eps_dyn1)-1 ;

save('K_series.mat', 'K_eps_dyn1', 'K_eps_dyn2', 'K_eps_csl')


figure(3) % IRFs for models
plot(0:T,K_eps_dyn1,0:T,K_eps_dyn2,'-.',0:T,K_eps_csl,'k--','LineWidth',2)
%axis([0 T+11 0.9*min([min(Kpath_TPI) min(Kpath_AMF)]) ...
%     1.1*max([max(Kpath_TPI) max(Kpath_AMF)])])

xlabel('period t','FontSize',15)
ylabel('aggregate capital K','FontSize',15)
title('IRFs for four solutions methods','FontSize',20)
legend('Dynare 1st order','Dynare 2nd order','CSL','Location','East')