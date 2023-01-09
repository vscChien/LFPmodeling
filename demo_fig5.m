% script for plotting figure 5
ref = load_data();
p   = load_sol();

%-----simulation-----
sim=cell(4,1); 
for s=1:4 % sites
  [sim{s}.MUA,sim{s}.CSD]=model(p(s,:),ref{s});
end

%-----plot MUA-----
figure('name','target MUA');
for s=1:4 % sites
   MUA{1}=ref{s}.MUA(201:400,:);
   MUA{2}=ref{s}.MUA(401:600,:);
   MUA{3}=ref{s}.MUA(1:200,:);
   MUA{4}=ref{s}.MUA(601:800,:);
   MUA{5}=ref{s}.MUA(801:1000,:);
   for i=1:5 % tones
       subplot(4,5,(s-1)*5+i);
       imagesc(1:200,1:16,MUA{i}');
       clim([-1 1]);colormap(flipud(coolwarm));
       hold on;
       plot(-MUA{i}*1.5+repmat(1:16,200,1),'k');
   end
end
set(gcf,'position',[0 0 800 600])


%-----plot CSD-----
figure('name','target CSD');
for s=1:4 % sites
   CSD{1}=ref{s}.CSD(201:400,:);
   CSD{2}=ref{s}.CSD(401:600,:);
   CSD{3}=ref{s}.CSD(1:200,:);
   CSD{4}=ref{s}.CSD(601:800,:);
   CSD{5}=ref{s}.CSD(801:1000,:);
   for i=1:5 % tones
       subplot(4,5,(s-1)*5+i);
       imagesc(1:200,3:14,CSD{i}');
       clim([-1 1]);colormap(flipud(coolwarm));
       hold on;
       plot(-CSD{i}*0.5+repmat(3:14,200,1),'k');
   end
end
set(gcf,'position',[0 0 800 600])



%-----plot simMUA-----
figure('name','simMUA');
for s=1:4 % sites
   simMUA=cell(5,1); % 5 tones
   simMUA{1}=sim{s}.MUA(201:400,:);
   simMUA{2}=sim{s}.MUA(401:600,:);
   simMUA{3}=sim{s}.MUA(1:200,:);
   simMUA{4}=sim{s}.MUA(601:800,:);
   simMUA{5}=sim{s}.MUA(801:1000,:);
   for i=1:5 % tones
       subplot(4,5,(s-1)*5+i);
       imagesc(1:200,1:16,simMUA{i}');
       clim([-1 1]);colormap(flipud(coolwarm));
       hold on;       
       plot(-simMUA{i}*1.5+repmat(1:16,200,1),'k');
   end
end
set(gcf,'position',[0 0 800 600])


%-----plot simCSD-----
figure('name','simCSD');
for s=1:4 % sites
   simCSD=cell(5,1); % 5 tones  
   simCSD{1}=sim{s}.CSD(201:400,:);
   simCSD{2}=sim{s}.CSD(401:600,:);
   simCSD{3}=sim{s}.CSD(1:200,:);
   simCSD{4}=sim{s}.CSD(601:800,:);
   simCSD{5}=sim{s}.CSD(801:1000,:);
   for i=1:5 % tones
       subplot(4,5,(s-1)*5+i);
       imagesc(1:200,3:14,simCSD{i}');
       clim([-1 1]);colormap(flipud(coolwarm));
       hold on;
       plot(-simCSD{i}*0.5+repmat(3:14,200,1),'k');
   end
end
set(gcf,'position',[0 0 800 600])



