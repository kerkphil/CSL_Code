% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% SIM_OUT.M produces output for SIMUL.M.  It is controlled 
% by several options to be described
% below.  All these options are set in OPTIONS.M,
% see that file for further details.
%
% The series will be plotted
% if SIM_GRAPH is set to 1.  In that case, further
% modifications can be done in particular with:
%   DO_ENLARGE : = 1, if you want large font sizes for the text on your plots.
%                     Good for slides.
%   SIM_JOINT  : = 1, if you want all series on the same graph, else = 0.
%   PRINT_FIG  : = 1, if you want plots to be printed on your printer
%   SAVE_FIG   : = 1, if you want plots to be saved as encaps. postscript. 
%                     Set PRINT_FIG = 0 also. The filenames are sim_ser1.eps, ...
%                     if SIM_JOINT = 0, and sim_data.eps is SIM_JOINT = 1.
%   SIM_PLOT_RAW : = 1, if you want a plot of the raw, unfiltered series, 
%                       even though you have chosen to HP-filter, DO_HP_FILTER = 1.                       
%                       Note, that if you have chosen to save figures, then
%                       the previously saved figures will be overwritten.
%                       This option is useful, if you want to look at the plot
%                       of the raw simulations, after having already seen the filtered
%                       ones: simply type
%                       SIM_PLOT_RAW = 1;
%                       sim_out;
%
% For printing the numbers of the autocorrelation table, 
% the following options are helpful:
%   
%  SIM_DO_DISP1: Set to = 1 to see printout of the autocorrelation matrix. 
%  SIM_DO_DISP2: Set to = 1 to see printout of the variance-covariance matrix.
%  SIM_DO_DISP3: Set to = 1 to see printout of the vector of variances.


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

n_select = max(size(SIM_SELECT));
if SIM_GRAPH,
   time_axis = (0:(SIM_LENGTH-1))/PERIOD + SIM_DATE0;
   if SIM_JOINT
      if SIM_PLOT_RAW,
         hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
            sim_raw(SIM_SELECT,1:min(SIM_LENGTH,SIM_MAX)));
         [mx,pos]=max(abs(sim_raw(SIM_SELECT,1:floor(SIM_CUT*min(SIM_LENGTH,SIM_MAX)))'));
         for var_index = 1:n_select,
%            text(time_axis(SIM_TXT_MARKER), sim_raw(var_index,SIM_TXT_MARKER),...
%              VARNAMES(SIM_SELECT(var_index),:));
            text(time_axis(pos(var_index)), sim_raw(var_index,pos(var_index)),...
              VARNAMES(SIM_SELECT(var_index),:));
         end;
      else
         hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
            sim_xyz(SIM_SELECT,1:min(SIM_LENGTH,SIM_MAX)));
         [mx,pos]=max(abs(sim_xyz(SIM_SELECT,1:floor(SIM_CUT*min(SIM_LENGTH,SIM_MAX)))'));
         for var_index = 1:n_select,
%            text(time_axis(SIM_TXT_MARKER), sim_xyz(var_index,SIM_TXT_MARKER),...
%              VARNAMES(SIM_SELECT(var_index),:));
            text(time_axis(pos(var_index)), sim_xyz(var_index,pos(var_index)),...
              VARNAMES(SIM_SELECT(var_index),:));
         end;
      end;
      set(hndl,'LineWidth',2);
      grid;
      if DO_HP_FILTER & ~SIM_PLOT_RAW,
         title('Simulated data (HP-filtered)');
      else
         title('Simulated data');
      end;
      xlabel('Year');
      ylabel('Percent deviation from steady state');
      if DO_ENLARGE,
         enlarge;
      end;
      if PRINT_FIG,
         disp('SIM_OUT: Printing plot of Simulated data');
         print;
      elseif SAVE_FIG
         disp('SIM_OUT: Saving plot of Simulated data.  Filename is sim_data.eps ...');
         print -deps sim_data.eps
      else                     
         disp('Inspect figure. Hit key when ready...');
         pause;
      end;
   else
      for var_index = 1:n_select,
         if SIM_PLOT_RAW,
            hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
               sim_raw(SIM_SELECT(var_index),1:min(SIM_LENGTH,SIM_MAX)));
            % text(time_axis(SIM_TXT_MARKER), sim_raw(var_index,SIM_TXT_MARKER),...
            %    VARNAMES(SIM_SELECT(var_index),:));
         else
            hndl = plot(time_axis(1:min(SIM_LENGTH,SIM_MAX)),...
               sim_xyz(SIM_SELECT(var_index),1:min(SIM_LENGTH,SIM_MAX)));
            % text(time_axis(SIM_TXT_MARKER), sim_xyz(var_index,SIM_TXT_MARKER),...
            %    VARNAMES(SIM_SELECT(var_index),:));
         end;
         set(hndl,'LineWidth',2);
         grid;
         if DO_HP_FILTER & ~SIM_PLOT_RAW,
            title(['Simulated data (HP-filtered): ',VARNAMES(SIM_SELECT(var_index),:)]);
         else
            title(['Simulated data: ',VARNAMES(SIM_SELECT(var_index),:)]);
         end;
         xlabel('Year');
         ylabel('Percent deviation from steady state');
         if DO_ENLARGE,
            enlarge;
         end;
         if PRINT_FIG,
            disp(['SIM_OUT.M: Printing simulation of ',VARNAMES(var_index,:),'...']);
            print;
         elseif SAVE_FIG
            if var_index == 1,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser1.eps ...'); 
               print -deps sim_ser1.eps
            elseif var_index == 2,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser2.eps ...'); 
               print -deps sim_ser2.eps
            elseif var_index == 3,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser3.eps ...'); 
               print -deps sim_ser3.eps
            elseif var_index == 4,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser4.eps ...'); 
               print -deps sim_ser4.eps
            elseif var_index == 5,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser5.eps ...'); 
               print -deps sim_ser5.eps
            elseif var_index == 6,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser6.eps ...'); 
               print -deps sim_ser6.eps
            elseif var_index == 7,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser7.eps ...'); 
               print -deps sim_ser7.eps
            elseif var_index == 8,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser8.eps ...'); 
               print -deps sim_ser8.eps
            elseif var_index == 9,
               disp(['SIM_OUT.M: Saving simulation of ',VARNAMES(var_index,:)]);
               disp( '         as encapsulated postscript file. Filename is sim_ser9.eps ...'); 
               print -deps sim_ser9.eps
            elseif var_index == 10,
               message = ['SIM_OUT.M: Sorry!  I can only save nine of your simulated series!      '
                          '           Perhaps be more selective about SIM_SELECT or change my code'];
               disp(message);
               warnings = [warnings;message];   
               disp(['           I cannot save the simulation of ',VARNAMES(var_index,:)]);
            elseif var_index > 10,
               disp(['           I cannot save the simulation of ',VARNAMES(var_index,:)]);
            end;
         else                     
            disp('Inspect figure. Hit key when ready...');
            pause;
         end;
      end;
   end;
end;
disp(' ');
if SIM_DO_DISP1 | SIM_DO_DISP2 | SIM_DO_DISP3,
   disp(        'SIM_OUT.M: Simulation-based calculation of moments');
   disp(sprintf('           Simulation length = %10d',SIM_LENGTH));
   if SIM_MODE == 2,
      disp(sprintf('           repeated %10d times.',SIM_N_SERIES));
   end;      
   disp('The variables are:');
   disp(VARNAMES(SIM_SELECT,:));
   disp(' ');
end;
if SIM_DO_DISP1,
   if DO_HP_FILTER,
      disp('Autocorrelation Table (HP-filtered series), corr(v(t+j),GNP(t)).  Last row shows j');
   else
      disp('Autocorrelation Table, corr(v(t+j),GNP(t)).  Last row shows j');
   end;
   for var_index = 1 : n_select,
      disp(sprintf('  %5.2f',autcor_sim(var_index,:)));
   end; 
   disp(sprintf('  %5.0f',autcor_sim(n_select+1,:)));
   disp(' ');
   if (SIM_N_SERIES > 3) & (SIM_MODE == 2),
      disp('Small sample standard errors for the Autocorrelation Table:');
      for var_index = 1 : n_select,
         disp(sprintf('  %5.2f',autcor_std(var_index,:)));
      end; 
      disp(sprintf('  %5.0f',autcor_sim(n_select+1,:)));
      disp(' ');
   end;
end;
if SIM_DO_DISP2,
   if DO_HP_FILTER,
      disp('Variance-Covariance Matrix (HP-filtered series):');
   else
      disp('Variance-Covariance Matrix:');
   end;
   for var_index = 1 : n_select,
      disp(sprintf(' %6.3f',covmat_sim(var_index,:)));
   end;
   disp(' ');
   if (SIM_N_SERIES > 3) & (SIM_MODE == 2),
      disp('Small sample standard errors for the Variance-Covariance Matrix:');
      for var_index = 1 : n_select,
         disp(sprintf(' %6.3f',covmat_std(var_index,:)));
      end;
      disp(' ');
   end;
end;
if SIM_DO_DISP3,
   if DO_HP_FILTER
      disp('Standard deviations (HP-filtered series):');
   else
      disp('Standard deviations:');
   end;
   disp(stdvec_sim);
   if (SIM_N_SERIES > 3) & (SIM_MODE == 2),
      disp('Small sample standard errors for the Standard deviations:');
      disp(stdvec_std);
   end;
end;
