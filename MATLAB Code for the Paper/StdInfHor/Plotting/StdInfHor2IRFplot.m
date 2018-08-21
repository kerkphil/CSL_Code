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

%-------------------------------------------------------------------------%
% Pull in 4 series
%-------------------------------------------------------------------------%
cd('../../StdInfHor/Dynare/1st order') ;
load('/Users/cc7768/Dropbox/Current State Linearization/MatLab/Std Inf Hor/Dynare/1st order/kmat2irf1ord.mat') ;
K_eps = kmat2
K_eps_dyn1 = [K_eps(1:end)]  ;
clear K_eps ;

cd('../../Dynare/2nd order') ;
load 'kmat2irf2ord.mat' ;
K_eps = kmat22
K_eps_dyn2 = [K_eps(1:end)]  ;
clear K_eps ;

cd('../../../StdInfHor/CSL') ;
load 'kmat2sdirfcsl.mat' Kcs ;
K_eps_csl = Kcs(1:end-1) ;
clear Kcs ;

cd('../../StdInfHor/Plotting') ;

T = length(K_eps_dyn1)-1 ;

save('K_2series.mat', 'K_eps_dyn1', 'K_eps_dyn2', 'K_eps_csl')

figure(4) % IRFs for models
plot(0:T,K_eps_dyn1,0:T,K_eps_dyn2,'-.',0:T,K_eps_csl,'k--','LineWidth',2)
%axis([0 T+11 0.9*min([min(Kpath_TPI) min(Kpath_AMF)]) ...
%     1.1*max([max(Kpath_TPI) max(Kpath_AMF)])])
xlabel('period t','FontSize',15)
ylabel('aggregate capital K','FontSize',15)
title('IRFs for four solutions methods','FontSize',20)
legend('Dynare 1st order','Dynare 2nd order','CSL','Location','East')