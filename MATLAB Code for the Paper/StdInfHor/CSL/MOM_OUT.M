% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% MOM_OUT produces output from the calculations done with MOMENTS.M,
% which is assumed to have been run just before.
% This program should be modified to suit tastes and needs.  Some options
% are given in the first few lines of this program.
% It is assumed that
% VARNAMES, a matrix with (m+n+k) rows, containing the variable names, has been set.
% The program overwrites freq1, diag_select, hndl, var_index


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

if DO_HP_GRAPH,
   freq1 = N_GRIDPOINTS/24;
   diag_select = (n_select+1)*(0:(n_select-1)) + 1; % Selects the diagonal
   hndl = plot(freqs(freq1:N_GRIDPOINTS/2), real(svv_raw(freq1:N_GRIDPOINTS/2,diag_select)));
   set(hndl,'LineWidth',2);
   for var_index = 1:n_select,
      text(freqs(MOM_TXT_MARKER), real(svv_raw(MOM_TXT_MARKER,diag_select(var_index))),...
           VARNAMES(HP_SELECT(var_index),:));
   end;
   title('Spectral densities, unfiltered');
   xlabel('Frequency');
   if DO_ENLARGE,
      enlarge;
   end;
   if PRINT_FIG,
      disp('MOM_OUT: Printing plot of unfiltered spectral densities');
      print;
   elseif SAVE_FIG
      disp('MOM_OUT: Saving plot of unfiltered spectral densities.  Filename is spec_raw.eps ...');
      print -deps spec_raw.eps
   else                     
      disp('Inspect figure. Hit key when ready...');
      pause;
   end;
   diag_select = (n_select+1)*(0:(n_select-1)) + 1; % Selects the diagonal
   hndl = plot(freqs(1:N_GRIDPOINTS/2), real(svv_fil(1:N_GRIDPOINTS/2,diag_select)));
   set(hndl,'LineWidth',2);
   for var_index = 1:n_select,
      text(freqs(MOM_TXT_MARKER), real(svv_fil(MOM_TXT_MARKER,diag_select(var_index))),...
           VARNAMES(HP_SELECT(var_index),:));
   end;
   title('Spectral densities, Hodrick-Prescott filtered');
   xlabel('Frequency');
   if DO_ENLARGE,
      enlarge;
   end;
   if PRINT_FIG,
      disp('MOM_OUT: Printing plot of filtered spectral densities');
      print;
   elseif SAVE_FIG
      disp('MOM_OUT: Saving plot of filtered spectral densities.  Filename is spec_fil.eps ...');
      print -deps spec_fil.eps
   else                     
      disp('Inspect figure. Hit key when ready...');
      pause;
   end;
end;
if DO_DISP1 | DO_DISP2 | DO_DISP3,
   disp(' ');
   disp('MOM_OUT.M: Frequency-domain method based calculation of moments');
   disp('The variables are:');
   disp(VARNAMES(HP_SELECT,:));
   disp(' ');
end;
if DO_HP_FILTER,
   if DO_DISP1,
      disp('Autocorrelation Table (HP-filtered series), corr(v(t+j),GNP(t)).  Last row shows j');
      for var_index = 1 : n_select,
         disp(sprintf('  %5.2f',autcor_fil(var_index,:)));
      end; 
      disp(sprintf('  %5.0f',autcor_fil(n_select+1,:)));
      disp(' ');
   end;
   if DO_DISP2,
      disp('Variance-Covariance Matrix, HP-filtered series:');
      for var_index = 1 : n_select,
         disp(sprintf(' %6.3f',covmat_fil(var_index,:)));
      end;
      disp(' ');
   end;
   if DO_DISP3,
      disp('Standard deviations, HP-filtered series:');
     disp(sqrt(varvec_fil));
   end;
else
   if DO_DISP1,
      disp('Autocorrelation Table (HP-filtered series), corr(v(t+j),GNP(t)).  Last row shows j');
      for var_index = 1 : n_select,
         disp(sprintf('  %5.2f',autcor_raw(var_index,:)));
      end; 
      disp(sprintf('  %5.0f',autcor_raw(n_select+1,:)));
      disp(' ');
   end;
   if DO_DISP2,
      disp('Variance-Covariance Matrix, HP-filtered series:');
      for var_index = 1 : n_select,
         disp(sprintf(' %6.3f',covmat_raw(var_index,:)));
      end;
      disp(' ');
   end;
   if DO_DISP3,
      disp('Standard deviations, HP-filtered series:');
     disp(sqrt(varvec_raw));
   end;
end;  
