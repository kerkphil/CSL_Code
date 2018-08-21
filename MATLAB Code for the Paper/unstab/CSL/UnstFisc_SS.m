function out = UnstFisc_SS(kbar,Bbar,zbar,param)

UnstFisc_dyn(in,param)

in = [kbar; Bbar; kbar; Bbar; kbar; Bbar; zbar; zbar];

temp = UnstFisc_dyn(in,param);

out = temp(1);