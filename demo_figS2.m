% script for plotting figure S2

s=3; % site 3
tmp=load(sprintf('.\\solution\\solution_nnmf_site%d.mat',s));

plotOn=1;
decomp_nnmf(tmp.p,tmp.ref,plotOn);


