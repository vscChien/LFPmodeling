% script for plotting figure 7a
ref = load_data();
p   = load_sol();

%-----plot simECD-----
figure('name','simECD');

yscale=[0.06,0.05,0.02,0.02];
for s=1:4 % sites
   [~,~,simECD,R]=model(p(s,:),ref{s});
   subp=[2,3,1,4,5];
   for i=1:5 % tones
       subplot(4,5,(s-1)*5+i);hold on;
       plot(R{subp(i)}.flow(:,1)*R{subp(i)}.dipoleInfo(1,3)+...
            R{subp(i)}.flow(:,4)*R{subp(i)}.dipoleInfo(4,3)+...
            R{subp(i)}.flow(:,8)*R{subp(i)}.dipoleInfo(8,3),'linewidth',2); % E1+E2+E3
       plot(R{subp(i)}.flow(:,2)*R{subp(i)}.dipoleInfo(2,3)+...
            R{subp(i)}.flow(:,5)*R{subp(i)}.dipoleInfo(5,3),'linewidth',2); % PV1+PV2
       plot(R{subp(i)}.flow(:,3)*R{subp(i)}.dipoleInfo(3,3)+...
            R{subp(i)}.flow(:,6)*R{subp(i)}.dipoleInfo(6,3),'linewidth',2); % PV1+PV2
       plot(R{subp(i)}.flow(:,7)*R{subp(i)}.dipoleInfo(7,3),'linewidth',2); % Th
       plot(R{subp(i)}.simECD,'k','linewidth',2);
       ylim([-1 1]*yscale(s))
       box on;
   end
end
set(gcf,'position',[0 0 800 600])