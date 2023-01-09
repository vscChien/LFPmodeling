% script for plotting figure 6
ref = load_data();
p   = load_sol();

%-----plot firing rates-----
figure('name','firing rates (column 1)');
for s=1:4 % sites
   [~,~,simMEG,Rall,g]=model(p(s,:),ref{s});  % simulation
   inputw=[p(s,25),p(s,26),1,p(s,27),p(s,28)];
   subp=[2,3,1,4,5];
   for i=1:5 % tones
       subplot(4,5,(s-1)*5+i);hold on;
       plot(Rall{subp(i)}.gS(1,:)+3,'linewidth',2,'color',[205,111,159]/256); % SOM1
       plot(Rall{subp(i)}.gP(1,:)+3,'linewidth',2,'color',[130,204,214]/256); % PV1
       plot(Rall{subp(i)}.gE(1,:)+3,'linewidth',2,'color',[1 1 1]*0.3);       % E1
       plot(Rall{subp(i)}.gE(3,:)+2,'linewidth',2,'color',[1 1 1]*0.3);       % E3
       plot(Rall{subp(i)}.gS(2,:)+1,'linewidth',2,'color',[205,111,159]/256); % SOM2  
       plot(Rall{subp(i)}.gP(2,:)+1,'linewidth',2,'color',[130,204,214]/256); % PV2
       plot(Rall{subp(i)}.gE(2,:)+1,'linewidth',2,'color',[1 1 1]*0.3);       % E2
       plot(Rall{subp(i)}.input*inputw(i),'linewidth',2,'color',[1 1 1]*0.3); % input
       plot(xlim,[1;1]*[1,2,3],'k')
       set(gca,'ytick',[0 1 2 3],'yticklabel',{'Th','L56','L4','L23'})
       ylim([0 4])
       box on;
   end
end
set(gcf,'position',[0 0 800 800])

