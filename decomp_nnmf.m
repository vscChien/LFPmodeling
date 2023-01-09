% Decompose MUA by non-negative matrix factorization (nnmf) 
%
% input:
%    p: free parameters [1x17]
%       p(1:8)   % synaptic kernel time constants
%       p(9:12)  % thalamic input strengths
%       p(13:17) % thalamic input decay levels 
%
%    ref: target data of the recording site
%       ref.CSD:     CSD   [200*5 timepoints x 12 channels]
%       ref.MUA:     MUA   [200*5 timepoints x 16 channels]
%       ref.times:   times [1 x 200 timepoints]
%       ref.freqs:   5 tones (BF tone and 4 non-BF tones) 
%       ref.delayIn: input delay (msec)
%       ref.MUApeak: peak value of MUA 
%
%    plotOn:  0: default
%             1: plot simulation regults
% output: 
%    simMUA: [200*5 timepoints x 16 channels] 
%    simCSD: [200*5 timepoints x 12 channels]
%    simECD: [200*5 timepoints x 1]
%    R:      simulation output
%
function [simMUA,simCSD,simECD,R]=decomp_nnmf(p,ref,plotOn)

    if nargin < 3
        plotOn=0;
    end
    p = rectify_p(p);

    [R,simMUA]=get_simMUA(ref);
    [R,simCSD]=get_simCSD(R,ref,p);
    [R,simECD]=get_simECD(R);
 
    %-----plot results-----
    if plotOn
        plot_mua_csd(R,simMUA,simCSD);
        plot_simECD(R,simECD);
        plot_rates_currents(R)
    end

end
%%=========================================================================
function p=rectify_p(p)
    p=p(:);
    boundary=zeros(17,2); % 17 parameters
    boundary(1:8,1)=1;    % min, synaptic kernel time constants
    boundary(1:8,2)=50;   % max, synaptic kernel time constants
    boundary(9:12,1)=0.1; % min, thalamic input strengths
    boundary(9:12,2)=1.2; % max, thalamic input strengths
    boundary(13:17,1)=0.1;% min, thalamic input decay levels
    boundary(13:17,2)=0.4;% max, thalamic input decay levels

    p(p<boundary(:,1))=boundary(p<boundary(:,1),1);
    p(p>boundary(:,2))=boundary(p>boundary(:,2),2); 
end
%%=========================================================================
function input=gen_input(ref,peak,ratio)        
   tauin=10;% msec
   input=ones(1,200); 
   input=ratio*input+...
         (1-ratio)*input.*exp(-(0:199)/tauin);
   
   input=[zeros(1,ref.delayIn),input];      
   input=input(1:200)*peak;  
   input=input';
end
%%=========================================================================
function [R,simMUA]=get_simMUA(ref)
    %-----decompose MUA to firing rates (nnmf)-----
    R.k=8; % number of components
    R.MUA=ref.MUA([201:600,1:200,601:1000],:); 
    R.CSD=ref.CSD([201:600,1:200,601:1000],:); 
    if isfield(ref,'rate')
       R.rate=ref.rate([201:600,1:200,601:1000],:);
       R.MUAsp=ref.MUAsp;
    else
      [R.rate,R.MUAsp] = nnmf(R.MUA,R.k-1);
    end
    simMUA=R.rate*R.MUAsp;
    R.R2mua=1-sumsqr(simMUA-R.MUA)/sumsqr(R.MUA-mean(R.MUA(:))); 
end

%%=========================================================================
function plot_mua_csd(R,simMUA,simCSD)
    figure('name','goodness of fit');
    subplot(121)
    h1=plot(R.MUA+repmat(16:-1:1,size(R.MUA,1),1),'r','linewidth',2);hold on
    h2=plot(simMUA+repmat(16:-1:1,size(R.MUA,1),1),'k','linewidth',2);
    ylim([0 17]);title(sprintf('MUA\n(%dcomps, R2=%g)',R.k,R.R2mua))
    set(gca,'ytick',1:16,'yticklabel',16:-1:1)
    plot([1;1]*(1:5)*200,ylim,'k')
    legend([h1(1),h2(1)],{'ref','sim'})
    xlabel('time (ms)');ylabel('channel')
    
    subplot(122)
    h1=plot(R.CSD+repmat(24:-2:2,size(R.CSD,1),1),'r','linewidth',2);hold on
    h2=plot(simCSD+repmat(24:-2:2,size(R.CSD,1),1),'k','linewidth',2);
    ylim([0 26]);title(sprintf('CSD\n(R2=%g)',R.R2csd))
    set(gca,'ytick',2:2:24,'yticklabel',14:-1:3)
    plot([1;1]*(1:5)*200,ylim,'k')
    legend([h1(1),h2(1)],{'ref','sim'})
    xlabel('time (ms)');
    set(gcf,'position',[0 0 800 800])
end
%%=========================================================================
function plot_simECD(R,simECD)
    figure('name','simECD');
    for i=2:R.k
        hold on;
        h1=plot(R.current(:,i)*R.sDipole(i),'color',[1 1 1]*0.7,'linewidth',2);   
    end
    h2=plot(R.current(:,1)*R.sDipole(1),'r','linewidth',2);
    h3=plot(simECD,'k','linewidth',2);
    title('simECD')
    xlabel('time (msec)');
    box on
    set(gcf,'position',[0 0 600 200])
    plot(xlim,[0;0],'k')
    plot([1;1]*(1:5)*200,ylim,'k')
    legend([h1,h2,h3],'nnmf','Th','sum')  
end
%%=========================================================================
function plot_rates_currents(R)
    areas=sum(R.MUAsp,2);
    figure('name','Details');
    k=R.k;
    for i=1:k
        subplot(k,4,i*4-3)
        plot(R.rate2(:,i),'k','linewidth',2);hold on;
        plot([1;1]*(1:5)*200,ylim,'k')
        title(sprintf('rate%d',i))
        if i==k,xlabel('time (ms)');end
        
        if i~=1
            subplot(k,4,i*4-2)
            bar(R.MUAsp(i-1,:),'edgecolor','k','facecolor','k')   
            set(gca,'xtick',[1,16])
            view([90 90]);
            title(sprintf('MUAsp%d (%g)',i,round(areas(i-1),2)))
        end
        if i==k,xlabel('channel');ylabel('weight');end
        
        subplot(k,4,i*4-1)
        plot(R.current(:,i),'k','linewidth',2);hold on;
        plot([1;1]*(1:5)*200,ylim,'k')
        title(sprintf('flow%d (%g ms)',i,R.tau(i)*1000))
        if i==k,xlabel('time (ms)');end
        
        subplot(k,4,i*4)
        bar(3:14,R.CSDsp(i,:),'k','facecolor','k')    
        xlim([2.5 14.5]);ylim([-1,1]*max(abs(R.CSDsp(i,:))));hold on;
        plot([R.cSink(i) R.cSource(i)]+2,[0 0],'r','linewidth',2)
        scatter(R.cSource(i)+2,0,'r_')
        set(gca,'xtick',[3,14])
        view([90 90])
        title(sprintf('CSDsp%d',i))
        if i==k,xlabel('channel');ylabel('weight');end
    end
    set(gcf,'position',[0 0 800 800])
end
%%=========================================================================
function [R,simCSD]=get_simCSD(R,ref,p)
    k=size(R.rate,2)+1;
    %----- gen input-----
    peak2=p(k+1);
    peak3=p(k+2);
    peak4=p(k+3);
    peak5=p(k+4);
    ratio1=p(k+5);
    ratio2=p(k+6);
    ratio3=p(k+7);
    ratio4=p(k+8);
    ratio5=p(k+9);

    input_BF     = gen_input(ref,1,ratio1);
    input_nonBF1 = gen_input(ref,peak2,ratio2);
    input_nonBF2 = gen_input(ref,peak3,ratio3);
    input_nonBF3 = gen_input(ref,peak4,ratio4);
    input_nonBF4 = gen_input(ref,peak5,ratio5);
    input=[input_nonBF1;input_nonBF2;input_BF;input_nonBF3;input_nonBF4];

    %------ synaptic kernels---------
    t  = (0:149)*1e-3;
    H  = exp(1)*1e-3; 
    tau= p(1:k)*1e-3;
    kernel=zeros(k,length(t));
    for i=1:k
      kernel(i,:)=H./tau(i)*t.*exp(-t/tau(i));
    end
    %-----currents-----
    R.current=zeros(size(R.rate,1),k);
    R.rate2=[input,R.rate];
    for i=1:k
      for tt=1:size(R.rate2,1)/200
        tmp = conv(R.rate2(tt*200-199:tt*200,i)',kernel(i,:)); 
        R.current(tt*200-199:tt*200,i)=tmp(1:200);
      end
    end
    R.tau=tau;
    R.CSDsp=pinv(R.current)*R.CSD;
    simCSD=R.current*R.CSDsp;
    R.R2csd=1-sumsqr(simCSD-R.CSD)/sumsqr(R.CSD-mean(R.CSD(:))); 
end
%%=========================================================================
function [R,simECD]=get_simECD(R)  
    nCh=size(R.CSDsp,2);
    tmp=R.CSDsp'; 
    posCSD=tmp.*(tmp>0);
    negCSD=-tmp.*(tmp<0);
    posCSD=posCSD./repmat(sum(posCSD,1),nCh,1);
    negCSD=negCSD./repmat(sum(negCSD,1),nCh,1);  
    R.cSink =   (1:nCh)*negCSD; % center of sink
    R.cSource = (1:nCh)*posCSD; % center of source
    R.sDipole = R.cSink-R.cSource;    % strength and direction of dipole
    R.cDipole = (R.cSink+R.cSource)/2; % center of dipole
    simECD=R.current*R.sDipole';
end
