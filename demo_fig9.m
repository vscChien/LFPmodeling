% script for plotting figure 9
ref = load_data();
p   = load_sol();

%-----plot parameters-----
figure('name','parameters');
for s=1:4 % sites
    
   freqs=ref{s}.freqs([2 3 1 4 5]); 
   lateralInhibition=p(s,[21 22 20 23 24]);
   subplot(2,1,1);semilogx(freqs,lateralInhibition,'-o','linewidth',2);hold on;
   grid on
   title('Lateral inhibition')
  
   inputStrength=[p(s,[25 26]), 1, p(s,[27 28])];
   subplot(2,1,2);semilogx(freqs,inputStrength,'-o','linewidth',2);hold on;
   grid on
   title('Input strength')
   xlabel('tone frequency (Hz)') 
end
set(gcf,'position',[0 0 800 300])


