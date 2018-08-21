%-------------------------------------------------------------------------%
% Compute solutions and IRFs to growth model with varying marginal tax
% rates
%-------------------------------------------------------------------------%
% This program calls the following .m, .mod, and .mat file(s)
%     * varyingMtax_par.m: Loads the parameters across models into memory
%     * varyingMtax_vfi.m: Value function iteration solution file
%     * varyingMtax_dy1.m: Dynare file computing 1st-order approximation of
%          policy functions
%     * varyingMtax_dy2.m: Dynare file computing 2nd-order approximation of
%          policy functions
%     * varyingMtax_csl.m: uses output from ClosedForm_d1.mod to make IRF
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/CSL/MatLab/varyingMtax') ;

%-------------------------------------------------------------------------%
% Run different solution methods
%-------------------------------------------------------------------------%

% Generate parameter values
run('../varyingMtax/ParamSet/varyingMtax_par.m') ;

% Run value function iteration solution
run('../varyingMtax/VFI/varyingMtax_vfi.m') ;

% Run Dynare 1st-order approximation
cd('Dynare1') ;
run('varyingMtax_dy1.m');
cd('../../varyingMtax/') ;

% Run Dynare 2nd-order approximation
cd('Dynare2') ;
run('varyingMtax_dy2.m')
cd('../../varyingMtax/') ;

% Run CSL approximation
cd('CSL') ;
run('varyingMtax_csl.m')
cd('../../varyingMtax/') ;

%-------------------------------------------------------------------------%
% Plot impulse response functions
%-------------------------------------------------------------------------%
clear ;
load '../varyingMtax/VFI/varyingMtax_vfi.mat' Tirf Tsim kirfvec_vfi1 ...
    kirfvec_vfi2 ksimvec_vfi3 eulerr_rms_vfi1 eulerr_rms_vfi2 ...
    eulerr_rms_vfi3 eulerr_mab_vfi1 eulerr_mab_vfi2 eulerr_mab_vfi3 ;
load '../varyingMtax/Dynare1/varyingMtax_dy1.mat' kirfvec_dy11 ...
    kirfvec_dy12 ksimvec_dy13 eulerr_rms_dy11 eulerr_rms_dy12 ...
    eulerr_rms_dy13 eulerr_mab_dy11 eulerr_mab_dy12 eulerr_mab_dy13 ...
    vfierr_rms_dy11 vfierr_rms_dy12 vfierr_rms_dy13 ;
load '../varyingMtax/Dynare2/varyingMtax_dy2.mat' kirfvec_dy21 ...
    kirfvec_dy22 ksimvec_dy23 eulerr_rms_dy21 eulerr_rms_dy22 ...
    eulerr_rms_dy23 eulerr_mab_dy21 eulerr_mab_dy22 eulerr_mab_dy23 ...
    vfierr_rms_dy21 vfierr_rms_dy22 vfierr_rms_dy23 ;
load '../varyingMtax/CSL/varyingMtax_csl.mat' kirfvec_csl1 ...
    kirfvec_csl2 ksimvec_csl3 eulerr_rms_csl1 eulerr_rms_csl2 ...
    eulerr_rms_csl3 eulerr_mab_csl1 eulerr_mab_csl2 eulerr_mab_csl3 ...
    vfierr_rms_csl1 vfierr_rms_csl2 vfierr_rms_csl3 ;

figure(1)
plot((0:1:Tirf),kirfvec_vfi1,(0:1:Tirf),kirfvec_dy11,(0:1:Tirf),...
    kirfvec_dy21,(0:1:Tirf),kirfvec_csl1) ;
xlabel('Period (t)','FontSize', 12);
ylabel('Aggregate Capital (K)','FontSize', 12);
%title('Comparison of methods on (kss,sigma) IRF','FontSize',18);
legend('VFI','Dyn1','Dyn2','CSL') ;

varyingMtaxEulErrs1 = [eulerr_rms_vfi1, eulerr_rms_dy11, eulerr_rms_dy21,...
                     eulerr_rms_csl1; ...
                     eulerr_mab_vfi1, eulerr_mab_dy11, eulerr_mab_dy21,...
                     eulerr_mab_csl1; ...
                     0, vfierr_rms_dy11, vfierr_rms_dy21, vfierr_rms_csl1] ;
display(varyingMtaxEulErrs1) ;

figure(2)
plot((0:1:Tirf),kirfvec_vfi2,(0:1:Tirf),kirfvec_dy12,(0:1:Tirf),...
    kirfvec_dy22,(0:1:Tirf),kirfvec_csl2);
xlabel('Period (t)','FontSize', 12);
ylabel('Aggregate Capital (K)','FontSize', 12);
%title('Comparison of methods on (1.1*kss,2*sigma) IRF','FontSize',18);
legend('VFI','Dyn1','Dyn2','CSL') ;

varyingMtaxEulErrs2 = [eulerr_rms_vfi2, eulerr_rms_dy12, eulerr_rms_dy22,...
                     eulerr_rms_csl2; ...
                     eulerr_mab_vfi2, eulerr_mab_dy12, eulerr_mab_dy22,...
                     eulerr_mab_csl2; ...
                     0, vfierr_rms_dy12, vfierr_rms_dy22, vfierr_rms_csl2] ;
display(varyingMtaxEulErrs2) ;

figure(3)
plot((0:1:Tsim),ksimvec_vfi3,(0:1:Tsim),ksimvec_dy13,(0:1:Tsim),...
    ksimvec_dy23,(0:1:Tsim),ksimvec_csl3);
xlabel('Period (t)','FontSize', 12);
ylabel('Aggregate Capital (K)','FontSize', 12);
%title('Comparison of methods on simulation','FontSize',18);
legend('VFI','Dyn1','Dyn2','CSL') ;

varyingMtaxEulErrs3 = [eulerr_rms_vfi3, eulerr_rms_dy13, eulerr_rms_dy23,...
                     eulerr_rms_csl3; ...
                     eulerr_mab_vfi3, eulerr_mab_dy13, eulerr_mab_dy23,...
                     eulerr_mab_csl3; ...
                     0, vfierr_rms_dy13, vfierr_rms_dy23, vfierr_rms_csl3] ;
display(varyingMtaxEulErrs3) ;

%-------------------------------------------------------------------------%
save varyingMtax_master.mat ;