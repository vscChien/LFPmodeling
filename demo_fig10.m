% script for plotting figure 10
ref = load_data();
p   = load_sol();

%-----calculate SSE and save to .mat file-----
% The scanning takes a few minutes.
for s=1:4
   filename=sprintf('.\\solution\\sensitivitycheck_site%d.mat',s);
   if ~exist(filename,'file')    
       fprintf('%s not found.\n',filename,s);
       p=[p0(s,:),1]; 
       % -----set scanning range-----
       boundary=zeros(29,2); % 28+1 parameters
       boundary(1:8,1)=0;      % min, connection strengths
       boundary(1:8,2)=15;     % max, connection strengths
       boundary([9 10],1)=0;    % min, ratio of thalamic inputs
       boundary([9 10],2)=10;   % max, ratio of thalamic inputs 
       boundary([11 12],1)=0;   % min, short-term plasticity 
       boundary([11 12],2)=4;   % max, short-term plasticity 
       boundary(13,1)=0.5;      % min, synaptic kernel time constants
       boundary(13,2)=1.5;      % max, synaptic kernel time constants
       boundary(14,1)=0.5;      % min, sigmoid function slopes
       boundary(14,2)=1.5;      % max, sigmoid function slopes
       boundary([15 16 17 18 19],1)=0;     % min, thalamic input decay levels
       boundary([15 16 17 18 19],2)=0.6;   % max, thalamic input decay levels
       boundary([20 21 22 23 24],1)=0;     % min, lateral inhibition E2→SOM1,2
       boundary([20 21 22 23 24],2)=15;    % max, lateral inhibition E2→SOM1,2
       boundary([25 26 27 28 29],1)=0.1;   % min, thalamic input strengths
       boundary([25 26 27 28 29],2)=1.2;   % max, thalamic input strengths
      
       target=[ref{s}.CSD(:);ref{s}.MUA(:)];
       SSE=cell(29,1);  % to be saved
       scanP=cell(29,1);% to be saved
       for idxP=1:29
           resolution=25;% scan resolution  
           scanP{idxP}=sort([linspace(boundary(idxP,1),boundary(idxP,2),resolution),p(idxP)]);
           tmpP=p;
           SSE{idxP}=zeros(length(scanP{idxP}),1);
           for i=1:length(scanP{idxP})       
               tmpP(idxP)=scanP{idxP}(i);      
               [simMUA,simCSD,~,~]=model(tmpP,ref{s});
               data=[simCSD(:);simMUA(:)];       
               SSE{idxP}(i)=sumsqr(target(:)-data(:));
               fprintf('site[%d/4]:SSE{%d/29}(%d/%d):%g\n',s,idxP,i,length(scanP{idxP}),SSE{idxP}(i));
           end
       end
       save(filename,'ref','p','boundary','scanP','SSE')
       fprintf('%s saved.\n',filename);
   end
end

%-----plot parameter sensitivity-----
figure('name','parameter sensitivity');
for s=1:4
   filename=sprintf('.\\solution\\sensitivitycheck_site%d.mat',s);
   clear ref p boundary scanP SSE
   load(filename);
   isubplot=[1:14,21,19:20,22:23, 27, 25:26,28:29,31:32, 34:35,33];
   for idxP=1:29
       nSSE=SSE{idxP};
       if 1
         nSSE=nSSE-min(nSSE);
         nSSE=nSSE/max(nSSE);
       end
       subplot(6,6,isubplot(idxP));hold on
       plot(scanP{idxP},nSSE,'-');
   end
end   
ttext={'Wee','Wpe','Wse','Wep','Wpp','Wsp','Wes','Wps',...
      'inE','inPV','adaptation','facilitation','tau','slope',...
      'decay1','decay2','decay3','decay4','decay5'...
      'lateral1','lateral2','lateral3','lateral4','lateral5'...
      'in2','in3','in4','in5','in1'};
for idxP=1:29
   subplot(6,6,isubplot(idxP));
   set(gca,'ColorOrderIndex',1);
   title(ttext{idxP});
end
for s=1:4
   filename=sprintf('.\\solution\\sensitivitycheck_site%d.mat',s); 
   clear ref p boundary scanP SSE
   load(filename);
   isubplot=[1:14,21,19:20,22:23, 27, 25:26,28:29,31:32, 34:35,33];
   for idxP=1:29
       subplot(6,6,isubplot(idxP));hold on     
       [~,ia,~]=intersect(scanP{idxP},p(idxP));
       h(s)=scatter(p(idxP),0,'^','filled');
       box on;
   end
end

set(gcf,'position',[0 0 1200 1000])

