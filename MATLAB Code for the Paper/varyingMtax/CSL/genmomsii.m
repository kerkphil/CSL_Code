function moms = genmomsii(Y,c,lagval)
%This function generates the moments that will be calculated and compared
%for CSL2012.  It will be saved in a file.

%Calculate AC, Var, StdErr for Y
ACY = corrcoef(Y(lagval + 1:end),Y(1:end- lagval));
VARY = var(Y);
STDERRY = sqrt(VARY);

%Calculate AC, Var, StdErr for c
ACC = corrcoef(c(2:end),c(1:end-1));
VARC = var(c);
STDERRC = sqrt(VARC);

moms = zeros(2,2);

moms(1,1) = ACY(1,2);
moms(2,1) = ACC(1,2);
moms(1,2) = STDERRY;
moms(2,2) = STDERRC;
% printmat(moms, 'Moments', 'Y C', 'AC STDERR')
%printmat(data, title, Rows, Columns)