% This file runs CSLClosedForm.m for each of the possible inputs and 
% saves the data and graphs in the folder K_irfs and the eul errs in
% the folder Eul_err_data

[eulerrs1, K_irf1] = CSLClosedForm('1stddev') ;
[eulerrs2, K_irf2] = CSLClosedForm('2stddev') ;
[eulerrsr, K_irfr] = CSLClosedForm('randsim') ;

save('Eul_err_data/eulerrs.mat', 'eulerrs1', 'eulerrs2', 'eulerrsr')

save('K_irfs/K_irf_data.mat', 'K_irf1', 'K_irf2', 'K_irfr')
