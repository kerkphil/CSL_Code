%-------------------------------------------------------------------------%
% Compute solutions and IRFs to standard infinite horizon problem for which
% the closed form solution is known
%-------------------------------------------------------------------------%
% This program calls the following .m, .mod, and .mat file(s)
%     * ClosedFormParams.m: Loads the parameters across models into memory
%     * ClosedForm_cs.m: Generates IRF from the closed form solution
%     * ClosedForm_d1.mod: Dynare file computing 1st-order approximation of
%          policy functions
%     * ClosedForm_d1run.m: uses output from ClosedForm_d1.mod to make IRF
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/ClosedForm') ;

%-------------------------------------------------------------------------%
% Run different solution methods
%-------------------------------------------------------------------------%

% Generate parameter values
% run('../ClosedForm/ParamSet/ClosedFormParams.m') ;

% Run closed form analytical solution
% run('../ClosedForm/ClosedFormSol/ClosedForm_ans.m') ;

% Run closed form value function iteration solution
% run('../ClosedForm/VFI/ClosedForm_vfi.m') ;

% Run Dynare 1st-order approximation of closed form solution
% cd('../ClosedForm/Dynare1/') ;
% clear ;
% ClosedForm_dy1 ;
% cd('../../ClosedForm/') ;

% Run Dynare 2nd-order approximation of closed form solution
% cd('../ClosedForm/Dynare2/') ;
% clear ;
% ClosedForm_dy2 ;
% cd('../../ClosedForm/') ;

% Run CSL approximation of closed form solution
% run('../ClosedForm/ClosedFormCSL/ClosedForm_csl.m') ;

%-------------------------------------------------------------------------%
% Plot impulse response functions
%-------------------------------------------------------------------------%
clear ;
load '../ClosedForm/ClosedFormSol/ClosedForm_ans.mat' Tirf Tsim ...
    kirfvec_ans1 kirfvec_ans2 ksimvec_ans3 eulerr_rms_ans1 ...
    eulerr_rms_ans2 eulerr_rms_ans3 eulerr_mab_ans1 eulerr_mab_ans2 ...
    eulerr_mab_ans3 ;
load '../ClosedForm/VFI/ClosedForm_vfi.mat' kirfvec_vfi1 kirfvec_vfi2 ...
    ksimvec_vfi3 eulerr_rms_vfi1 eulerr_rms_vfi2 eulerr_rms_vfi3 ...
    eulerr_mab_vfi1 eulerr_mab_vfi2 eulerr_mab_vfi3 anlerr_rms_vfi1 ...
    anlerr_rms_vfi2 anlerr_rms_vfi3 ;
load '../ClosedForm/Dynare1/ClosedForm_dy1.mat' kirfvec_dy11 ...
    kirfvec_dy12 ksimvec_dy13 eulerr_rms_dy11 eulerr_rms_dy12 ...
    eulerr_rms_dy13 eulerr_mab_dy11 eulerr_mab_dy12 eulerr_mab_dy13 ...
    anlerr_rms_dy11 anlerr_rms_dy12 anlerr_rms_dy13 ;
load '../ClosedForm/Dynare2/ClosedForm_dy2.mat' kirfvec_dy21 ...
    kirfvec_dy22 ksimvec_dy23 eulerr_rms_dy21 eulerr_rms_dy22 ...
    eulerr_rms_dy23 eulerr_mab_dy21 eulerr_mab_dy22 eulerr_mab_dy23 ...
    anlerr_rms_dy21 anlerr_rms_dy22 anlerr_rms_dy23 ;
load '../ClosedForm/ClosedFormCSL/ClosedForm_csl.mat' kirfvec_csl1 ...
    kirfvec_csl2 ksimvec_csl3 eulerr_rms_csl1 eulerr_rms_csl2 ...
    eulerr_rms_csl3 eulerr_mab_csl1 eulerr_mab_csl2 eulerr_mab_csl3 ...
    anlerr_rms_csl1 anlerr_rms_csl2 anlerr_rms_csl3 ;

figure(1)
plot((0:1:Tirf),kirfvec_ans1,(0:1:Tirf),kirfvec_vfi1,(0:1:Tirf),kirfvec_dy11,(0:1:Tirf),kirfvec_dy21,(0:1:Tirf),kirfvec_csl1)
xlabel('Period (t)','FontSize', 12);
ylabel('Aggregate Capital (K)','FontSize', 12);
%title('Comparison of methods on (kss,sigma) IRF','FontSize',18);
legend('Analytical Sol.','VFI','Dyn1','Dyn2','CSL') ;

ClosedFormEulErrs1 = [eulerr_rms_ans1, eulerr_rms_vfi1, eulerr_rms_dy11, eulerr_rms_dy21, eulerr_rms_csl1; ...
                      eulerr_mab_ans1, eulerr_mab_vfi1, eulerr_mab_dy11, eulerr_mab_dy21, eulerr_mab_csl1; ...
                      0, anlerr_rms_vfi1, anlerr_rms_dy11, anlerr_rms_dy21, anlerr_rms_csl1] ;
display(ClosedFormEulErrs1) ;

figure(2)
plot((0:1:Tirf),kirfvec_ans2,(0:1:Tirf),kirfvec_vfi2,(0:1:Tirf),kirfvec_dy12,(0:1:Tirf),kirfvec_dy22,(0:1:Tirf),kirfvec_csl2)
xlabel('Period (t)','FontSize', 12);
ylabel('Aggregate Capital (K)','FontSize', 12);
%title('Comparison of methods on (1.1*kss,2*sigma) IRF','FontSize',18);
legend('Analytical Sol.','VFI','Dyn1','Dyn2','CSL') ;

ClosedFormEulErrs2 = [eulerr_rms_ans2, eulerr_rms_vfi2, eulerr_rms_dy12, eulerr_rms_dy22, eulerr_rms_csl2; ...
                      eulerr_mab_ans2, eulerr_mab_vfi2, eulerr_mab_dy12, eulerr_mab_dy22, eulerr_mab_csl2; ...
                      0, anlerr_rms_vfi2, anlerr_rms_dy12, anlerr_rms_dy22, anlerr_rms_csl2] ;
display(ClosedFormEulErrs2) ;

ClosedFormEulErrs3 = [eulerr_rms_ans3, eulerr_rms_vfi3, eulerr_rms_dy13, eulerr_rms_dy23, eulerr_rms_csl3; ...
                      eulerr_mab_ans3, eulerr_mab_vfi3, eulerr_mab_dy13, eulerr_mab_dy23, eulerr_mab_csl3; ...
                      0, anlerr_rms_vfi3, anlerr_rms_dy13, anlerr_rms_dy23, anlerr_rms_csl3] ;
display(ClosedFormEulErrs3) ;

%-------------------------------------------------------------------------%
save ClosedFormMaster.mat ;
