% script for plotting figure 7b

%-----run nnmf-----
simECD=cell(4,1);
current=cell(4,1);
sDipole=cell(4,1);
for s=1:4 % sites
    tmp=load(sprintf('.\\solution\\solution_nnmf_site%d.mat',s)); 
    [~,~,simECD{s},R]=decomp_nnmf(tmp.p,tmp.ref,0);
    current{s}=R.current;
    sDipole{s}=R.sDipole;
end

%-----plot simECD by nnmf-----
figure('name','simECD by nnmf');

yscale=[0.08,0.06 0.05 0.02];
for s=1:4 % sites
   for i=1:5 % tones
       subplot(4,5,(s-1)*5+i);hold on;
       for c=2:8
         h=plot(current{s}(1+200*(i-1):200*i,c)*sDipole{s}(c),'linewidth',2,'color',[1 1 1]*0.7);
       end
       h1=plot(current{s}(1+200*(i-1):200*i,1)*sDipole{s}(1),'r','linewidth',2);
       h2=plot(simECD{s}(1+200*(i-1):200*i,1),'k','linewidth',2);
       ylim([-1 1]*yscale(s))
       box on;
   end
end
legend([h(1),h1,h2],'nnmf','Th','sum')
set(gcf,'position',[0 0 800 600])

