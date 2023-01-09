% script for plotting figure S1
ref = load_data();
p   = load_sol();

%--------------------------------------
for s=1:4 % sites
    [simMUA,simCSD]=model(p(s,:),ref{s}); % simulation 
    %-----plot MUA&CSD-----  
    figname=sprintf('site %d',s);
    figure('name',figname);
    timeindex=[201:600,1:200,601:1000];
    subplot(121);hold on
    h1=plot(ref{s}.MUA(timeindex,:)+repmat(16:-1:1,1000,1),'r','linewidth',2);
    h2=plot(simMUA(timeindex,:)+repmat(16:-1:1,1000,1),'k','linewidth',2);
    set(gca,'ytick',1:16,'yticklabel',16:-1:1,'xtick',200,'xticklabel',200);
    ylim([0 17]);
    plot([1;1]*[200,400,600,800],ylim,'k');
    plot(xlim,[1;1]*(1:16),'k')
    ylabel('channel');xlabel('time (ms)')
    R2=1-sumsqr(ref{s}.MUA(:)-simMUA(:))/sumsqr(ref{s}.MUA(:)-mean(ref{s}.MUA(:))); 
    title(sprintf('MUA (R^2=%g)',R2));
    legend([h1(1),h2(1)],'ref','sim')
    box on;
    %--------------
    subplot(122);hold on
    h1=plot(ref{s}.CSD(timeindex,:)/2+repmat(14:-1:3,1000,1),'r','linewidth',2);
    h2=plot(simCSD(timeindex,:)/2+repmat(14:-1:3,1000,1),'k','linewidth',2);
    set(gca,'ytick',3:14,'yticklabel',14:-1:3,'xtick',200,'xticklabel',200);
    ylim([2 15]);
    plot([1;1]*[200,400,600,800],ylim,'k');
    plot(xlim,[1;1]*(3:14),'k')
    ylabel('channel');xlabel('time (ms)')
    R2=1-sumsqr(ref{s}.CSD(:)-simCSD(:))/sumsqr(ref{s}.CSD(:)-mean(ref{s}.CSD(:))); 
    title(sprintf('CSD (R^2=%g)',R2));
    legend([h1(1),h2(1)],'ref','sim')
    box on;
    set(gcf,'position',[0 0 800 800])
end


