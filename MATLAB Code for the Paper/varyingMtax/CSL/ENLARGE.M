% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% ENLARGE.M modifies an existing plot with a graphics handle
% called hndl by enlarging the character sizes.  This is useful
% for projecting the graphs via an overhead projector, for example.
% The procedure is called as part of impresp.m
% Options like font sizes are set in OPTIONS.M
% 
% Note: the size of the letters are controlled by LET_SIZE for titles
% etc. and LAB_SIZE for labels.  If you want to change them BEFORE you have
% run enlarge.m for a particular graph, just change LET_SIZE and LAB_SIZE
% as you please: you can do this most easily by modifying options.m.
% The default setting by MATLAB is 12.
% If you want to change them AFTER you
% have run enlarge once already (say, after looking at a graph and not
% liking it), you need to type in the following lines:
%
% OLD_LET_SIZE = LET_SIZE;
% OLD_LAB_SIZE = LAB_SIZE;
% LET_SIZE = <whatever new size you prefer>;
% LAB_SIZE = <whatever new size you prefer>;
% enlarge;
%
% IF that didn't work for some reason, try the following:
%
% par = get(hndl(1),'Parent');
% tit = get(par,'Title');
% OLD_LET_SIZE = get(tit,'FontSize');
% chil = get(par,'Children');
% OLD_LAB_SIZE = get(chil(1),'FontSize');
% LET_SIZE = <whatever new size you prefer>;
% LAB_SIZE = <whatever new size you prefer>;
% enlarge;
%
% To make enlarge work for the next graph you look at, you may then have to reset:
%
% OLD_LET_SIZE = 12;
% OLD_LAB_SIZE = 12;

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.


% Options are set in OPTIONS.M
par = get(hndl(1),'Parent');
txt_hndl = findobj(par,'FontSize',OLD_LET_SIZE);
chil = get(par,'Children');
label_chil = findobj(chil,'FontSize',OLD_LAB_SIZE);
for j = 1 : max(size(txt_hndl)),
   set(txt_hndl(j),'FontSize',LET_SIZE);
end;
for j = 1 : max(size(label_chil)),
   set(label_chil(j),'FontSize',LAB_SIZE);
end;
