% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% IMPRESP.M calculates and plots impulse responses to shocks one percent
% in size to each one of the exogenous variables z(t), keeping all
% others at zero in the first period.  It is assumed,
% that SOLVE.M has been executed before, so that the matrices
% NN, PP, QQ, RR and SS are available, describing the law of motion
%   x(t) = PP x(t-1) + QQ z(t)
%   y(t) = RR x(t-1) + SS z(t)
%   z(t) = NN z(t-1) + epsilon(t)
% The following options should have been chosen beforehand (see also file OPTIONS.M)
%   HORIZON    : how far out should the impulse responses be calculated
%   PERIOD     : number of periods per year, i.e. 12 for monthly, 4 for quarterly
%   DO_PLOTS   : = 1, if plots should be made, = 0, if not.
%   IMP_SELECT : a vector containing the indices of the variables to be
%             plotted.  To plot all variables, set IMP_SELECT = 1 : (m+n+k),
%             where m=dim(x), n=dim(y) and k=dim(z);
%   VARNAMES   : an array with (m+n+k) rows, containing the variable names.
%   The following options are set in OPTIONS.M:
%   TXT_MARKER : a number indicating where the label for the responses should be written.
%   DO_ENLARGE : = 1, if you want large font sizes for the text on your plots.  Good for slides.
%   PRINT_FIG  : = 1, if you want plots to be printed on your printer
%   SAVE_FIG   : = 1, if you want plots to be saved as encaps. postscript.  Set PRINT_FIG = 0 also.
% The program calculates:
% Response: the response of all variables x(t), y(t), z(t) to each shock.
%   Response is of size (m_states+n_endog+k_exog)*HORIZON.
%   Since Response is overwritten, each time a new shock is analyzed, 
%   the results are collected in 
% Resp_mat = [ Response to first shock
%              Response to second shock
%              ...                     ]
%
% The program also defines 
% Time_axis,m_states,n_endog,k_exog,II_contemp,II_lag,hndl,
% thus overwriting variables with these names that might have been used before.


% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

% Most options are now set in OPTIONS.M

% Calculations
[m_states,k_exog] = size(QQ);
[n_endog,k_exog]  = size(SS);
Time_axis = (0 : (HORIZON-1))/PERIOD;
Resp_mat = [];
for shock_counter = 1 : k_exog,
   Response = zeros(m_states+n_endog+k_exog,HORIZON);
   Response(m_states+n_endog+shock_counter,1) = 1;
   II_lag = [ PP, zeros(m_states,n_endog),zeros(m_states,k_exog)
              RR, zeros(n_endog, n_endog),zeros(n_endog, k_exog)
              zeros(k_exog,(m_states+n_endog)), NN                ];
   II_contemp = eye(m_states+n_endog+k_exog) + ...
        [ zeros(m_states,(m_states+n_endog)), QQ
          zeros(n_endog, (m_states+n_endog)), SS
          zeros(k_exog,  (m_states+n_endog)), zeros(k_exog,k_exog) ];
   % describing [x(t)',y(t)',z(t)']'= II_contemp*II_lag*[x(t-1)',y(t-1)',z(t-1)']';
   Response(:,1) = II_contemp*Response(:,1);
   for time_counter = 2 : HORIZON,
      Response(:,time_counter) = II_contemp*II_lag*Response(:,time_counter-1);
   end;
   Resp_mat = [ Resp_mat 
                Response ];
   if DO_PLOTS,
      hndl = plot(Time_axis,0*Time_axis, ...
                  Time_axis,Response(IMP_SELECT,:));
      set(hndl(2:max(size(hndl))),'LineWidth',2);
      grid;
      title(['Impulse responses to a shock in ', ...
          VARNAMES(m_states+n_endog+shock_counter,:)]);
      xlabel('Years after shock');
      ylabel('Percent deviation from steady state');
      for varindex = IMP_SELECT,
         text(Time_axis(min([TXT_MARKER,HORIZON/2])), ...
            Response(varindex,min([TXT_MARKER,HORIZON/2])),...
            VARNAMES(varindex,:));
      end;
      if DO_ENLARGE,
         enlarge;
      end;
      if PRINT_FIG,
         disp(['IMPRESP.M: Printing response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:),'...']);
         print;
      elseif SAVE_FIG
         if shock_counter == 1,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp1.eps ...'); 
            print -deps impresp1.eps
         elseif shock_counter == 2,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp2.eps ...'); 
            print -deps impresp2.eps
         elseif shock_counter == 3,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp3.eps ...'); 
            print -deps impresp3.eps
         elseif shock_counter == 4,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp4.eps ...'); 
            print -deps impresp4.eps
         elseif shock_counter == 5,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp5.eps ...'); 
            print -deps impresp5.eps
         elseif shock_counter == 6,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp6.eps ...'); 
            print -deps impresp6.eps
         elseif shock_counter == 7,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp7.eps ...'); 
            print -deps impresp7.eps
         elseif shock_counter == 8,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp8.eps ...'); 
            print -deps impresp8.eps
         elseif shock_counter == 9,
            disp(['IMPRESP.M: Saving response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:)]);
            disp( '         as encapsulated postscript file. Filename is impresp9.eps ...'); 
            print -deps impresp9.eps
         elseif shock_counter > 9,
            disp('IMPRESP.M: Sorry!  I can only save nine of your impulse responses!');
            disp(['         So, response to a shock in ',VARNAMES(m_states+n_endog+shock_counter,:),' is not saved.']);
            disp('         If you want me to save more, you need to change my code.');
         end;
      else                     
         disp('Inspect figure. Hit key when ready...');
         pause;
      end;
   end;
end;
