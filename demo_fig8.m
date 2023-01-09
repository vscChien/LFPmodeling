% script for plotting figure 8
ref = load_data();
p   = load_sol();

%-----simulation-----
R=cell(4,1);
for s=1:4
   [~,~,~,R{s},g]=model(p(s,:),ref{s});
end

%-----plot MUA spatial profile-----
figure('name','MUAsp');
ttext={'E1','E2','E3','PV1','PV2','SOM1','SOM2'};
for s=1:4 % sites
       Rtmp=R{s}{1};
       for i=1:7 % populations
           subplot(7,4,(i-1)*4+s)  
           bar(Rtmp.MUAsp(i,:),'k');
           xlim([0.5 16.5]);ylim([0,1]);
           set(gca,'xtick',[1,16])
           view([90 90])
           title(ttext{i})
       end
end
set(gcf,'position',[0 0 400 800])


%-----plot CSD spatial profile-----
figure('name','CSDsp');
ttext={'E1','P1','S1','E2','P2','S2','Th','E3'};
for s=1:4 % sites
       Rtmp=R{s}{1};
       subpl=[2,5,7,3,6,8,1,4];
       for i=1:8 % source inputs
           subplot(8,4,(subpl(i)-1)*4+s)                      
           bar(3:14,Rtmp.CSDsp(:,i)/max(abs(Rtmp.CSDsp(:))),'k');
           hold on;
           plot([Rtmp.dipoleInfo(i,1);Rtmp.dipoleInfo(i,2)]+2,[0;0],'r','linewidth',2);
           scatter(Rtmp.dipoleInfo(i,1)+2,0,'r_');
           xlim([2.5 14.5]);ylim([-1,1]);
           set(gca,'xtick',[3,14])
           view([90 90])
           title(ttext{i});
       end
end
set(gcf,'position',[0 0 400 800])
