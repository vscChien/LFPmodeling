% model: total 7 populations
%        E*3: pyramidal cell
%        S*2: SOM interneuron
%        P*2: PV interneuron
%
% input:
%    p: free parameters [1x28]
%       p(1:8)   % connection strengths 
%       p(9:10)  % ratio of thalamic inputs
%       p(11:12) % short-term plasticity 
%       p(13)    % synaptic kernel time constants
%       p(14)    % sigmoid function slopes
%       p(15:19) % thalamic input decay levels
%       p(20:24) % lateral inhibition E2→SOM1,2
%       p(25:28) % thalamic input strengths
%
%    ref: target data of the recording site
%       ref.CSD:     CSD   [200*5 timepoints x 12 channels]
%       ref.MUA:     MUA   [200*5 timepoints x 16 channels]
%       ref.times:   times [1 x 200 timepoints]
%       ref.freqs:   5 tones (BF tone and 4 non-BF tones) 
%       ref.delayIn: input delay (msec)
%       ref.MUApeak: peak value of MUA 
%
%    returnW: 0: default
%             1: return g without simulation
%
%    plotOn:  0: default
%             1: plot simulation regults
%
% output: 
%    simMUA: [200*5 timepoints x 16 channels] 
%    simCSD: [200*5 timepoints x 12 channels]
%    simECD: [200*5 timepoints x 1]
%    R:      [5 tones x 1] cells, simulation output
%    g:      configuration of the simulation
%
function [simMUA,simCSD,simECD,R,g]=model(p,ref,returnW,plotOn)
 
    if nargin < 3
        returnW=0;
    end
    if nargin < 4
        plotOn=0;
    end   
    if length(p)<29
        p(29)=1; % p(29) is only used in sensitivity check (demo_fig10.m) 
    end  
    p(p<0)=0; % p should be larger than 0.

    %-----configuration----- 
    g = model_setting(p,ref);
    if returnW==1 
       simMUA=[];simCSD=[];R=[];simECD=[];
       return
    end

    %-----run simulation-----
    R=cell(5,1); % 5 tones
    % BF condition
    g.ratio = p(15);  
    g       = gen_input(g);
    g.Wse   = g.Wse_bf;
    g.Wet   = g.Wet_bf;
    g.Wpt   = g.Wpt_bf;
    g.Wst   = g.Wst_bf;
    R{1}    = model_core(g);

    % non-BF1 condition
    g.ratio = p(16);  
    g       = gen_input(g);
    g.Wse   = g.Wse_nbf1; 
    g.Wet   = g.Wet_nbf1; 
    g.Wpt   = g.Wpt_nbf1; 
    g.Wst   = g.Wst_nbf1; 
    R{2}    = model_core(g);

    % non-BF2 condition
    g.ratio = p(17);
    g       = gen_input(g);
    g.Wse   = g.Wse_nbf2; 
    g.Wet   = g.Wet_nbf2; 
    g.Wpt   = g.Wpt_nbf2; 
    g.Wst   = g.Wst_nbf2; 
    R{3}    = model_core(g);
    
    % non-BF3 condition
    g.ratio = p(18); 
    g       = gen_input(g);
    g.Wse   = g.Wse_nbf3; 
    g.Wet   = g.Wet_nbf3;
    g.Wpt   = g.Wpt_nbf3; 
    g.Wst   = g.Wst_nbf3; 
    R{4}    = model_core(g);
    
    % non-BF4 condition
    g.ratio = p(19); 
    g       = gen_input(g);
    g.Wse   = g.Wse_nbf4; 
    g.Wet   = g.Wet_nbf4; 
    g.Wpt   = g.Wpt_nbf4;
    g.Wst   = g.Wst_nbf4; 
    R{5}    = model_core(g); 
    
    %-----simMUA, simCSD, and simECD-----
    [R,simMUA] = get_simMUA(R,ref);
    [R,simCSD] = get_simCSD(R,ref);    
    [R,simECD] = get_simECD(R);
   
    %-----plot results-----
    if plotOn
        plot_rates_column1(R,p);
        plot_rates_column2(R);
        plot_mua_csd(ref,simMUA,simCSD);
        plot_mua_csd_detail(R{1});
    end
end
%%========================================================================= 
function g=gen_input(g)        
 
   input=ones(g.nTone,g.soundLength); 
   input=g.ratio*input+...
         (1-g.ratio)*input.*exp(-(0:g.soundLength-1)/g.tauin);
   
   input=[zeros(g.nTone,g.onset+g.delay),input];      
   g.input=input(:,1:g.soundLength);
   g.dt=1e-4;
   g.times=(1:size(g.input,2))*g.dt;
    
end

%%=========================================================================
function R = model_core(g)

    %-----connection strengths-----
    Wee=g.Wee;
    Wse=g.Wse;
    Wpe=g.Wpe;
   
    Wep=g.Wep;
    Wsp=g.Wsp;
    Wpp=g.Wpp;
   
    Wes=g.Wes;
    Wps=g.Wps;
  
    Wet=g.Wet;
    Wst=g.Wst;
    Wpt=g.Wpt;    
    
    %-----synaptic kernel time constants----
    tauee1=g.taue1(1); 
    taueeNMDA1=g.taueNMDA1;
    tause1=g.taue1(2);
    taupe1=g.taue1(3);

    tauep1=g.taup1(1); 
    tausp1=g.taup1(2);
    taupp1=g.taup1(3);

    taues1=g.taus1(1); 
    tauesGABAB1=g.tausGABAB1;
    taups1=g.taus1(2);
    
    tauee2=g.taue2(1); 
    taueeNMDA2=g.taueNMDA2;
    tause2=g.taue2(2);
    taupe2=g.taue2(3);

    tauep2=g.taup2(1); 
    tausp2=g.taup2(2);
    taupp2=g.taup2(3);

    taues2=g.taus2(1); 
    tauesGABAB2=g.tausGABAB2;
    taups2=g.taus2(2);
    
    %-----synaptic gains-----
    Hee=g.He(1); 
    HeeNMDA=g.HeNMDA;
    Hse=g.He(2);
    Hpe=g.He(3);

    Hep=g.Hp(1); 
    Hsp=g.Hp(2);
    Hpp=g.Hp(3);

    Hes=g.Hs(1); 
    HesGABAB=g.HsGABAB;
    Hps=g.Hs(2);
      
    %-----background---
    B=g.B; 
    
    %-----sigmoid-----
    theta_e=g.theta_e;
    slope_e=g.slope_e;  
    theta_s=g.theta_s;
    slope_s=g.slope_s; 
    theta_p=g.theta_p;
    slope_p=g.slope_p; 
    
    %-----STP-----
    STPee  = g.STPee; 
    STPse  = g.STPse; 
    STPth  = g.STPth; 

    Uee    = 1;      
    tauree = g.tauree;
    taufee = 1000000;
    kree   = g.kree;
    kfee   = 0;

    Use    = g.Use;  
    taurse = 1000000;
    taufse = g.taufse;
    krse   = 0;
    kfse   = g.kfse;
    
    Uth    = 1;      
    taurth = g.taurth;
    taufth = 1000000;
    krth   = g.krth;
    kfth   = 0;

    %-----input-----
    input  = g.input; 
    dt     = g.dt;
    soundLength=size(input,2); 

    %-----variables initialization----- 
    Eee=zeros(g.Ne,soundLength);      % EPSP from E (E population)
    EeeNMDA=zeros(g.Ne,soundLength);  % EPSP from E (E population)
    Ees=zeros(g.Ne,soundLength);      % IPSP from SOM (E population)
    EesGABAB=zeros(g.Ne,soundLength); % IPSP from SOM (E population)
    Eep=zeros(g.Ne,soundLength);      % IPSP from PV (E population)
    Sse=zeros(g.Ns,soundLength);      % EPSP from E  (SOM population)
    Ssp=zeros(g.Ns,soundLength);      % IPSP from PV (SOM population)
    Ppe=zeros(g.Np,soundLength);      % EPSP from E  (PV population)
    Pps=zeros(g.Np,soundLength);      % IPSP from SOM (PV population)
    Ppp=zeros(g.Np,soundLength);      % IPSP from PV (PV population)
    
    Eee2=zeros(g.Ne,soundLength);     % 2nd EPSP from E (E population)
    EeeNMDA2=zeros(g.Ne,soundLength); % 2nd EPSP from E (E population)
    Ees2=zeros(g.Ne,soundLength);     % 2nd IPSP from SOM (E population)
    EesGABAB2=zeros(g.Ne,soundLength);% 2nd IPSP from SOM (E population)
    Eep2=zeros(g.Ne,soundLength);     % 2nd IPSP from PV (E population)
    Sse2=zeros(g.Ns,soundLength);     % 2nd EPSP from E  (SOM population)
    Ssp2=zeros(g.Ns,soundLength);     % 2nd IPSP from PV (SOM population)
    Ppe2=zeros(g.Np,soundLength);     % 2nd EPSP from E  (PV population)
    Pps2=zeros(g.Np,soundLength);     % 2nd IPSP from SOM (PV population)
    Ppp2=zeros(g.Np,soundLength);     % 2nd IPSP from PV (PV population)  
    
    gE=zeros(g.Ne,soundLength);       % firing rate (E population)
    gS=zeros(g.Ns,soundLength);       % firing rate (SOM population)    
    gP=zeros(g.Np,soundLength);       % firing rate (PV population) 
    
    xee=ones(g.Ne,soundLength);       % STP; release probability
    xse=ones(g.Ne,soundLength);
    xth=ones(g.Nt,soundLength);
    uee=ones(g.Ne,soundLength);       % STP; utility
    use=ones(g.Ne,soundLength); 
    uth=ones(g.Nt,soundLength);

    if STPee, uee=uee*Uee;  end          
    if STPse, use=use*Use;  end
    if STPth, uth=uth*Uth;  end
    
    %-----specific EPSPs (stand-alone)-----
    EA=zeros(g.Ne,soundLength);  % EPSP from E1 on {E1;E2;E3}
    EA2=zeros(g.Ne,soundLength); % 2nd  EPSP  
    EB=zeros(g.Ne,soundLength);  % EPSP from E2 on {E1;E2;E3}
    EB2=zeros(g.Ne,soundLength); % 2nd EPSP   
    EC=zeros(g.Ne,soundLength);  % EPSP from E3 on {E1;E2;E3}
    EC2=zeros(g.Ne,soundLength); % 2nd EPSP 
    ED=zeros(g.Ne,soundLength);  % EPSP from Thamus on {E1;E2;E3}
    ED2=zeros(g.Ne,soundLength); % 2nd EPSP
    
    ENMDAA=zeros(g.Ne,soundLength);  % EPSP from E1 on {E1;E2;E3}
    ENMDAA2=zeros(g.Ne,soundLength); % 2nd  EPSP  
    ENMDAB=zeros(g.Ne,soundLength);  % EPSP from E2 on {E1;E2;E3}
    ENMDAB2=zeros(g.Ne,soundLength); % 2nd EPSP   
    ENMDAC=zeros(g.Ne,soundLength);  % EPSP from E3 on {E1;E2;E3}
    ENMDAC2=zeros(g.Ne,soundLength); % 2nd EPSP 
    ENMDAD=zeros(g.Ne,soundLength);  % EPSP from Thamus on {E1;E2;E3}
    ENMDAD2=zeros(g.Ne,soundLength); % 2nd EPSP
      
    EepA=zeros(g.Ne,soundLength);    % IPSP from PV1 on {E1;E2;E3}
    EepA2=zeros(g.Ne,soundLength);   % 2nd  IPSP  
    EepB=zeros(g.Ne,soundLength);    % IPSP from PV2 on {E1;E2;E3}
    EepB2=zeros(g.Ne,soundLength);   % 2nd IPSP
    
    EesA=zeros(g.Ne,soundLength);    % IPSP from SOM1 on {E1;E2;E3}
    EesA2=zeros(g.Ne,soundLength);   % 2nd  IPSP  
    EesB=zeros(g.Ne,soundLength);    % IPSP from SOM2 on {E1;E2;E3}
    EesB2=zeros(g.Ne,soundLength);   % 2nd IPSP
    
    EesGABABA=zeros(g.Ne,soundLength);  % IPSP from SOM1 on {E1;E2;E3}
    EesGABABA2=zeros(g.Ne,soundLength); % 2nd  IPSP  
    EesGABABB=zeros(g.Ne,soundLength);  % IPSP from SOM2 on {E1;E2;E3}
    EesGABABB2=zeros(g.Ne,soundLength); % 2nd IPSP
    
   %-----simulation-----
    for t=1:soundLength-1                
        %-----potential to rate----- 
        gE(:,t)=sigm(Eee(:,t)*g.AMPA_ratio + EeeNMDA(:,t)*(1-g.AMPA_ratio) - Ees(:,t)*(1-g.GABAB_ratio) -EesGABAB(:,t)*g.GABAB_ratio -Eep(:,t),theta_e,slope_e); 
        gS(:,t)=sigm(Sse(:,t)-Ssp(:,t),theta_s,slope_s); 
        gP(:,t)=sigm(Ppe(:,t)-Pps(:,t)-Ppp(:,t),theta_p,slope_p);
      
        %-----rate to potential----- 
        % EPSP from E (E population) AMPA
        dEee  = Eee2(:,t);
        dEee2 = Hee*(Wee*(xee(:,t).*uee(:,t).*gE(:,t)) + Wet*(xth(:,t).*uth(:,t).*input(:,t)) + B) - (tauee2+tauee1)*Eee2(:,t)/tauee2/tauee1 - Eee(:,t)/tauee2/tauee1;
        % EPSP from E (E population) NMDA 
        dEeeNMDA  = EeeNMDA2(:,t);
        dEeeNMDA2 = HeeNMDA*(Wee*(xee(:,t).*uee(:,t).*gE(:,t)) + Wet*(xth(:,t).*uth(:,t).*input(:,t)*g.Th_NMDA_ratio) + B) - (taueeNMDA2+taueeNMDA1)*EeeNMDA2(:,t)/taueeNMDA2/taueeNMDA1 - EeeNMDA(:,t)/taueeNMDA2/taueeNMDA1;                    
        % IPSP from SOM (E population) slow GABA_A
        dEes  = Ees2(:,t);
        dEes2 = Hes*Wes*gS(:,t) - (taues2+taues1)*Ees2(:,t)/taues2/taues1 - Ees(:,t)/taues2/taues1;
        % IPSP from SOM (E population) GABA_B
        dEesGABAB  = EesGABAB2(:,t);
        dEesGABAB2 = HesGABAB*Wes*gS(:,t) - (tauesGABAB2+tauesGABAB1)*EesGABAB2(:,t)/tauesGABAB2/tauesGABAB1 - EesGABAB(:,t)/tauesGABAB2/tauesGABAB1;           
        % IPSP from PV (E population)
        dEep  = Eep2(:,t);
        dEep2 = Hep*Wep*gP(:,t) - (tauep2+tauep1)*Eep2(:,t)/tauep2/tauep1 - Eep(:,t)/tauep2/tauep1;
        % EPSP from E  (SOM population)
        dSse  = Sse2(:,t); 
        dSse2 = Hse*(Wse*(xse(:,t).*use(:,t).*gE(:,t)) + Wst*(xth(:,t).*uth(:,t).*input(:,t))) - (tause2+tause1)*Sse2(:,t)/tause2/tause1 - Sse(:,t)/tause2/tause1;
        % IPSP from PV (SOM population)
        dSsp  = Ssp2(:,t);
        dSsp2 = Hsp*Wsp*gP(:,t) - (tausp2+tausp1)*Ssp2(:,t)/tausp2/tausp1 - Ssp(:,t)/tausp2/tausp1;
        % EPSP from E  (PV population)
        dPpe  = Ppe2(:,t);
        dPpe2 = Hpe*(Wpe*gE(:,t) + Wpt*(xth(:,t).*uth(:,t).*input(:,t))) - (taupe2+taupe1)*Ppe2(:,t)/taupe2/taupe1 - Ppe(:,t)/taupe2/taupe1;  
        % IPSP from SOM (PV population)
        dPps  = Pps2(:,t);
        dPps2 = Hps*Wps*gS(:,t) - (taups2+taups1)*Pps2(:,t)/taups2/taups1 - Pps(:,t)/taups2/taups1;
        % IPSP from PV (PV population)    
        dPpp  = Ppp2(:,t);
        dPpp2 = Hpp*Wpp*gP(:,t) - (taupp2+taupp1)*Ppp2(:,t)/taupp2/taupp1 - Ppp(:,t)/taupp2/taupp1;
        
        
        Eee(:,t+1)=Eee(:,t)+dEee*dt;
        EeeNMDA(:,t+1)=EeeNMDA(:,t)+dEeeNMDA*dt; 
        Ees(:,t+1)=Ees(:,t)+dEes*dt;
        EesGABAB(:,t+1)=EesGABAB(:,t)+dEesGABAB*dt;
        Eep(:,t+1)=Eep(:,t)+dEep*dt; 
        Sse(:,t+1)=Sse(:,t)+dSse*dt; 
        Ssp(:,t+1)=Ssp(:,t)+dSsp*dt; 
        Ppe(:,t+1)=Ppe(:,t)+dPpe*dt; 
        Pps(:,t+1)=Pps(:,t)+dPps*dt; 
        Ppp(:,t+1)=Ppp(:,t)+dPpp*dt; 
        

        Eee2(:,t+1)=Eee2(:,t)+dEee2*dt; 
        EeeNMDA2(:,t+1)=EeeNMDA2(:,t)+dEeeNMDA2*dt; 
        Ees2(:,t+1)=Ees2(:,t)+dEes2*dt; 
        EesGABAB2(:,t+1)=EesGABAB2(:,t)+dEesGABAB2*dt; 
        Eep2(:,t+1)=Eep2(:,t)+dEep2*dt; 
        Sse2(:,t+1)=Sse2(:,t)+dSse2*dt; 
        Ssp2(:,t+1)=Ssp2(:,t)+dSsp2*dt; 
        Ppe2(:,t+1)=Ppe2(:,t)+dPpe2*dt; 
        Pps2(:,t+1)=Pps2(:,t)+dPps2*dt; 
        Ppp2(:,t+1)=Ppp2(:,t)+dPpp2*dt; 
        
        %-----STP-----
        if STPee
           dxee=(1-xee(:,t))/tauree - kree*uee(:,t).*xee(:,t).*gE(:,t);
           duee=(Uee-uee(:,t))/taufee + kfee*Uee*(1-uee(:,t)).*gE(:,t);              
           xee(:,t+1)=xee(:,t)+dxee*dt;
           uee(:,t+1)=uee(:,t)+duee*dt; 
        end 
        if STPse
           dxie=(1-xse(:,t))/taurse - krse*use(:,t).*xse(:,t).*gE(:,t);
           duie=(Use-use(:,t))/taufse + kfse*Use*(1-use(:,t)).*gE(:,t);              
           xse(:,t+1)=xse(:,t)+dxie*dt;
           use(:,t+1)=use(:,t)+duie*dt; 
        end
        if STPth
           dxth=(1-xth(:,t))/taurth - krth*uth(:,t).*xth(:,t).*input(:,t);
           duth=(Uth-uth(:,t))/taufth + kfth*Uth*(1-uth(:,t)).*input(:,t);              
           xth(:,t+1)=xth(:,t)+dxth*dt;
           uth(:,t+1)=uth(:,t)+duth*dt; 
        end
        
        %-----------------specific EPSPs (stand-alone)-------------
        maskA=[1,0,0;... % E1 to E1 
               1,0,0;... %    to E2 
               1,0,0];   %    to E3
        maskB=[0,1,0;... % E2 to E1
               0,1,0;... %    to E2
               0,1,0];   %    to E3
        maskC=[0,0,1;... % E3 to E1 
               0,0,1;... %    to E2 
               0,0,1];   %    to E3
        maskA = kron([1,0;0,0],maskA); % two-column 
        maskB = kron([1,0;0,0],maskB); % two-column 
        maskC = kron([1,0;0,0],maskC); % two-column 
        
        dEA  = EA2(:,t);
        dEA2 = Hee*((maskA.*Wee)*(xee(:,t).*uee(:,t).*gE(:,t))) - (tauee2+tauee1)*EA2(:,t)/tauee2/tauee1 - EA(:,t)/tauee2/tauee1;   
        
        dEB  = EB2(:,t);
        dEB2 = Hee*((maskB.*Wee)*(xee(:,t).*uee(:,t).*gE(:,t))) - (tauee2+tauee1)*EB2(:,t)/tauee2/tauee1 - EB(:,t)/tauee2/tauee1;

        dEC  = EC2(:,t);
        dEC2 = Hee*((maskC.*Wee)*(xee(:,t).*uee(:,t).*gE(:,t))) - (tauee2+tauee1)*EC2(:,t)/tauee2/tauee1 - EC(:,t)/tauee2/tauee1;
             
        dED  = ED2(:,t);
        dED2 = Hee*(Wet*(xth(:,t).*uth(:,t).*input(:,t))) - (tauee2+tauee1)*ED2(:,t)/tauee2/tauee1 - ED(:,t)/tauee2/tauee1;


        EA(:,t+1)=EA(:,t)+dEA*dt;
        EA2(:,t+1)=EA2(:,t)+dEA2*dt;        
        EB(:,t+1)=EB(:,t)+dEB*dt;
        EB2(:,t+1)=EB2(:,t)+dEB2*dt;       
        EC(:,t+1)=EC(:,t)+dEC*dt;
        EC2(:,t+1)=EC2(:,t)+dEC2*dt;      
        ED(:,t+1)=ED(:,t)+dED*dt;
        ED2(:,t+1)=ED2(:,t)+dED2*dt;
 
        dENMDAA  = ENMDAA2(:,t);
        dENMDAA2 = HeeNMDA*((maskA.*Wee)*(xee(:,t).*uee(:,t).*gE(:,t))) - (taueeNMDA2+taueeNMDA1)*ENMDAA2(:,t)/taueeNMDA2/taueeNMDA1 - ENMDAA(:,t)/taueeNMDA2/taueeNMDA1;   
        
        dENMDAB  = ENMDAB2(:,t);
        dENMDAB2 = HeeNMDA*((maskB.*Wee)*(xee(:,t).*uee(:,t).*gE(:,t))) - (taueeNMDA2+taueeNMDA1)*ENMDAB2(:,t)/taueeNMDA2/taueeNMDA1 - ENMDAB(:,t)/taueeNMDA2/taueeNMDA1;

        dENMDAC  = ENMDAC2(:,t);
        dENMDAC2 = HeeNMDA*((maskC.*Wee)*(xee(:,t).*uee(:,t).*gE(:,t))) - (taueeNMDA2+taueeNMDA1)*ENMDAC2(:,t)/taueeNMDA2/taueeNMDA1 - ENMDAC(:,t)/taueeNMDA2/taueeNMDA1;
             
        dENMDAD  = ENMDAD2(:,t);
        dENMDAD2 = HeeNMDA*(Wet*(xth(:,t).*uth(:,t).*input(:,t))*g.Th_NMDA_ratio) - (taueeNMDA2+taueeNMDA1)*ENMDAD2(:,t)/taueeNMDA2/taueeNMDA1 - ENMDAD(:,t)/taueeNMDA2/taueeNMDA1;

        ENMDAA(:,t+1)=ENMDAA(:,t)+dENMDAA*dt;
        ENMDAA2(:,t+1)=ENMDAA2(:,t)+dENMDAA2*dt;        
        ENMDAB(:,t+1)=ENMDAB(:,t)+dENMDAB*dt;
        ENMDAB2(:,t+1)=ENMDAB2(:,t)+dENMDAB2*dt;       
        ENMDAC(:,t+1)=ENMDAC(:,t)+dENMDAC*dt;
        ENMDAC2(:,t+1)=ENMDAC2(:,t)+dENMDAC2*dt;      
        ENMDAD(:,t+1)=ENMDAD(:,t)+dENMDAD*dt;
        ENMDAD2(:,t+1)=ENMDAD2(:,t)+dENMDAD2*dt;

        
        % IPSP from PV1, PV2 (E population)
        maskepA=[1,0;... % PV1 to E1 
                 1,0;... %     to E2 
                 1,0];   %     to E3
        maskepB=[0,1;... % PV2 to E1
                 0,1;... %     to E2
                 0,1];   %     to E3
        maskepA = kron([1,0;0,0],maskepA); % two-column   
        maskepB = kron([1,0;0,0],maskepB); % two-column
             
        dEepA  = EepA2(:,t);
        dEepA2 = Hep*(maskepA.*Wep)*gP(:,t) - (tauep2+tauep1)*EepA2(:,t)/tauep2/tauep1 - EepA(:,t)/tauep2/tauep1;
        
        dEepB  = EepB2(:,t);
        dEepB2 = Hep*(maskepB.*Wep)*gP(:,t) - (tauep2+tauep1)*EepB2(:,t)/tauep2/tauep1 - EepB(:,t)/tauep2/tauep1;
              
        
        EepA(:,t+1)=EepA(:,t)+dEepA*dt;
        EepA2(:,t+1)=EepA2(:,t)+dEepA2*dt;      
        EepB(:,t+1)=EepB(:,t)+dEepB*dt;
        EepB2(:,t+1)=EepB2(:,t)+dEepB2*dt;
        
        % IPSP from SOM1, SOM2 (E population)
        maskesA=[1,0;... % SOM1 to E1 
                 1,0;... %      to E2 
                 1,0];   %      to E3
        maskesB=[0,1;... % SOM2 to E1
                 0,1;... %      to E2
                 0,1];   %      to E3
        maskesA = kron([1,0;0,0],maskesA); % two-column   
        maskesB = kron([1,0;0,0],maskesB); % two-column     
             
        dEesA  = EesA2(:,t);
        dEesA2 = Hes*(maskesA.*Wes)*gS(:,t) - (taues2+taues1)*EesA2(:,t)/taues2/taues1 - EesA(:,t)/taues2/taues1;
        
        dEesB  = EesB2(:,t);
        dEesB2 = Hes*(maskesB.*Wes)*gS(:,t) - (taues2+taues1)*EesB2(:,t)/taues2/taues1 - EesB(:,t)/taues2/taues1;
      
        EesA(:,t+1)=EesA(:,t)+dEesA*dt;
        EesA2(:,t+1)=EesA2(:,t)+dEesA2*dt;      
        EesB(:,t+1)=EesB(:,t)+dEesB*dt;
        EesB2(:,t+1)=EesB2(:,t)+dEesB2*dt;
        
        dEesGABABA  = EesGABABA2(:,t);
        dEesGABABA2 = HesGABAB*(maskesA.*Wes)*gS(:,t) - (tauesGABAB2+tauesGABAB1)*EesGABABA2(:,t)/tauesGABAB2/tauesGABAB1 - EesGABABA(:,t)/tauesGABAB2/tauesGABAB1;
        
        dEesGABABB  = EesGABABB2(:,t);
        dEesGABABB2 = HesGABAB*(maskesB.*Wes)*gS(:,t) - (tauesGABAB2+tauesGABAB1)*EesGABABB2(:,t)/tauesGABAB2/tauesGABAB1 - EesGABABB(:,t)/tauesGABAB2/tauesGABAB1;
      
        EesGABABA(:,t+1)=EesGABABA(:,t)+dEesGABABA*dt;
        EesGABABA2(:,t+1)=EesGABABA2(:,t)+dEesGABABA2*dt;      
        EesGABABB(:,t+1)=EesGABABB(:,t)+dEesGABABB*dt;
        EesGABABB2(:,t+1)=EesGABABB2(:,t)+dEesGABABB2*dt;        
    end
        
    %-----return output----- 
    R.Nt=g.Nt;
    R.Ne=g.Ne;
    R.Ns=g.Ns;
    R.Np=g.Np;
    R.Wee=Wee;
    R.Wse=Wse;
    R.Wpe=Wpe;
    R.Wep=Wep;
    R.Wsp=Wsp;
    R.Wpp=Wpp;
    R.Wes=Wes;
    R.Wps=Wps;
    R.Wet=Wet;
    R.Wst=Wst;
    R.Wpt=Wpt;
    
    R.Eee=Eee(:,1:10:end)*g.AMPA_ratio + EeeNMDA(:,1:10:end)*(1-g.AMPA_ratio); %{E1,E2,E3,Th} to {E1;E2;E3}
    R.EA =EA(:,1:10:end)*g.AMPA_ratio + ENMDAA(:,1:10:end)*(1-g.AMPA_ratio);   %{E1} to {E1;E2;E3}
    R.EB =EB(:,1:10:end)*g.AMPA_ratio + ENMDAB(:,1:10:end)*(1-g.AMPA_ratio);   %{E2} to {E1;E2;E3}
    R.EC =EC(:,1:10:end)*g.AMPA_ratio + ENMDAC(:,1:10:end)*(1-g.AMPA_ratio);   %{E3} to {E1;E2;E3}
    R.ED =ED(:,1:10:end)*g.AMPA_ratio + ENMDAD(:,1:10:end)*(1-g.AMPA_ratio);   %{Th} to {E1;E2;E3}
    R.EepA=EepA(:,1:10:end); % PV1 to {E1;E2;E3}
    R.EepB=EepB(:,1:10:end); % PV2 to {E1;E2;E3}
    R.EesA=EesA(:,1:10:end)*(1-g.GABAB_ratio) + EesGABABA(:,1:10:end)*g.GABAB_ratio; % SOM1 to {E1;E2;E3}
    R.EesB=EesB(:,1:10:end)*(1-g.GABAB_ratio) + EesGABABB(:,1:10:end)*g.GABAB_ratio; % SOM2 to {E1;E2;E3}
    R.Ees=Ees(:,1:10:end)*(1-g.GABAB_ratio) + EesGABAB(:,1:10:end)*g.GABAB_ratio;
    R.Eep=Eep(:,1:10:end);
    R.Sse=Sse(:,1:10:end);
    R.Ssp=Ssp(:,1:10:end);
    R.Ppe=Ppe(:,1:10:end);
    R.Pps=Pps(:,1:10:end);
    R.Ppp=Ppp(:,1:10:end);
       
    R.E=Eee*g.AMPA_ratio+EeeNMDA*(1-g.AMPA_ratio)...
        -Ees*(1-g.GABAB_ratio)-EesGABAB*g.GABAB_ratio...
        -Eep; R.E=R.E(:,1:10:end);
    R.S=Sse-Ssp;     R.S=R.S(:,1:10:end);
    R.P=Ppe-Pps-Ppp; R.P=R.P(:,1:10:end);
    R.gE=gE(:,1:10:end);
    R.gS=gS(:,1:10:end);
    R.gP=gP(:,1:10:end);
    R.cee=uee.*xee; R.cee=R.cee(:,1:10:end);
    R.cse=use.*xse; R.cse=R.cse(:,1:10:end);
    R.cth=uth.*xth; R.cth=R.cth(:,1:10:end);
    R.input=g.input(:,1:10:end);
    
    R.flow=[sum(R.EA(1:3,:),1);...  % 1: mean EPSP from E1 
            sum(R.EepA(1:3,:),1);...% 2: mean IPSP from PV1
            sum(R.EesA(1:3,:),1);...% 3: mean IPSP from SOM1 
            sum(R.EB(1:3,:),1);...  % 4: mean EPSP from E2        
            sum(R.EepB(1:3,:),1);...% 5: IPSP from PV2
            sum(R.EesB(1:3,:),1);...% 6: IPSP from SOM2
            sum(R.ED(1:3,:),1);...  % 7: mean EPSP from Th
            sum(R.EC(1:3,:),1)]';   % 8: mean EPSP from E3
    R.rates=[R.gE(1:3,:)',R.gP(1:2,:)',R.gS(1:2,:)'];     
     
end

%%=========================================================================
function [R,simMUA]=get_simMUA(R,ref) 

    rate=[R{1}.rates(:,1:7);R{2}.rates(:,1:7);R{3}.rates(:,1:7);...
          R{4}.rates(:,1:7);R{5}.rates(:,1:7)];

    % cellDensity*maxFiring
    % [ref] M. Beierleinl, J.R. Gibson, B.W. Connors (2003) 
    % [ref] Keller, D., Erö, C., & Markram, H. (2018)
    MUAspLim=[59.4*128400/3;... %E1
              59.4*128400/3;... %E2
              59.4*128400/3;... %E3
              271.7*4345/2;...  %PV1
              271.7*4345/2;...  %PV2
              120.7*2142/2;...  %SOM1
              120.7*2142/2];    %SOM2
    MUAspLim=MUAspLim/max(MUAspLim(:));   % constraint: sum=1   
    MUA=ref.MUA(:); %[200*16 x 1]
    nCh=16;
    nPop=7;
    option=  optimset('display','off');
    

    rate=kron(eye(nCh),rate);%[200*16 x 7*16]
    % constraints
    Aeq=kron(ones(1,nCh),eye(nPop));
    beq=MUAspLim;%[7x1]
    lb=zeros(nPop*nCh,1);
    ub=[];
    %-----gradient descent to find best scalar-----
    eta=1; % learning rate
    opsteps = 10;
    criteria=1e-3;
    scalar=zeros(opsteps,1);
    resnorm=zeros(opsteps,1);
        
    scalar(1)=1;% intitial
    h=1e-6;     % inital step
           
    for j = 1:opsteps         
        [x,resnorm(j),res1]=lsqlin(rate,MUA,[],[],Aeq,beq*scalar(j),lb,ub,[],option);
        [~,~,res2]=lsqlin(rate,MUA,[],[],Aeq,beq*(scalar(j)+h),lb,ub,[],option);
        
        J = (res2-res1)/h;
        d = -eta*pinv(J'*J)*J'*res1;  
        scalar(j+1) = max(scalar(j)+d,0);
      
        if (j>=2) && (abs(resnorm(j)-resnorm(j-1))< criteria)
            break;
        end
    end

    simMUA=reshape(rate*x,1000,nCh);      
    MUAsp=reshape(x,nPop,nCh);         

    for i=1:5
       R{i}.MUAsp=MUAsp;
    end
    R{1}.simMUA=simMUA(1:200,:);
    R{2}.simMUA=simMUA(201:400,:);
    R{3}.simMUA=simMUA(401:600,:);    
    R{4}.simMUA=simMUA(601:800,:);
    R{5}.simMUA=simMUA(801:1000,:);
end

%==========================================================================
function [Rall,simCSD]=get_simCSD(Rall,ref)

   flow=[Rall{1}.flow; Rall{2}.flow; Rall{3}.flow; Rall{4}.flow; Rall{5}.flow];
   try
       L=flow; 
       y0=ref.CSD;       
       x=mvregress(L,y0); 
       xPos=x; xPos(xPos<0)=0; % positive part
       xNeg=x; xNeg(xNeg>0)=0; % negative part
       xPos=xPos./repmat(sum(xPos,2),[1 size(xPos,2)]); xPos(isnan(xPos))=0;
       xNeg=-xNeg./repmat(sum(xNeg,2),[1 size(xNeg,2)]); xNeg(isnan(xNeg))=0;
       b2=xPos+xNeg;
       b2(isnan(b2))=0; % avoid devide-by-zero
   catch
       keyboard;
   end
   CSDsp0=xPos+xNeg;
   CSDsp0=CSDsp0/max(max(L*CSDsp0));
   CSDsp0(isnan(CSDsp0))=0;
  
   simCSD=L*CSDsp0;
   for i=1:5
     Rall{i}.CSDsp=b2';     
     Rall{i}.CSDsp0=CSDsp0';
   end
   
   Rall{1}.simCSD=simCSD(1:200,:);
   Rall{2}.simCSD=simCSD(201:400,:);  
   Rall{3}.simCSD=simCSD(401:600,:); 
   Rall{4}.simCSD=simCSD(601:800,:); 
   Rall{5}.simCSD=simCSD(801:1000,:);     
end

%%=========================================================================
function [R,simECD]=get_simECD(R)

   nCh=size(R{1}.CSDsp,1);
   tmp=R{1}.CSDsp; 
   posCSD=tmp.*(tmp>0);
   negCSD=-tmp.*(tmp<0);
 
   cSink = (1:nCh)*negCSD;   % center of sink
   cSource = (1:nCh)*posCSD; % center of source
   sDipole = cSink-cSource;    % strength and direction of dipole
   cDipole = (cSink+cSource)/2; % center of dipole
                           
   for i=1:length(R)  
      R{i}.dipoleInfo=[cSink',cSource',sDipole'];
      R{i}.simECD=R{i}.flow*sDipole';     
   end
   simECD=[R{1}.simECD;R{2}.simECD;R{3}.simECD;R{4}.simECD;R{5}.simECD];
end

%%=========================================================================
function m=sigm(v, theta, slope)
   m = max(1./(1+exp(slope.*(theta-v))) - 1/(1+exp(slope*(theta))),0); 
end

%%=========================================================================
function plot_rates_column1(R,p)
    inputs=[1,p(25),p(26),p(27),p(28)];
    figure;
    idx=[2,3,1,4,5];
    for i=1:5

        subplot(3,5,i)% rate(L23)
        hold on;
        plot(R{idx(i)}.gE(1,:),'k','linewidth',2);
        plot(R{idx(i)}.gP(1,:),'c','linewidth',2);
        plot(R{idx(i)}.gS(1,:),'r','linewidth',2);
        plot(R{idx(i)}.cee(1,:),'g');
        plot(R{idx(i)}.cse(1,:),'b');
        ylim([0 1])
        legend('E1','PV1','SOM1','dep','fac')
        title('firing rate(L23)')
        grid on;box on;

        subplot(3,5,i+5)% rate(L4)
        hold on;
        plot(R{idx(i)}.gE(3,:),'k','linewidth',2);
        plot(R{idx(i)}.cee(3,:),'g');
        plot(R{idx(i)}.cse(3,:),'b');
        plot(R{idx(i)}.input(1,:)*inputs(idx(i)),'k');
        ylim([0 1])
        legend('E3','dep','fac','input')
        title('firing rate(L4)')
        grid on;box on;

        subplot(3,5,i+10)% rate(L5)
        hold on;
        plot(R{idx(i)}.gE(2,:),'k','linewidth',2);
        plot(R{idx(i)}.gP(2,:),'c','linewidth',2);
        plot(R{idx(i)}.gS(2,:),'r','linewidth',2);
        plot(R{idx(i)}.cee(2,:),'g');
        plot(R{idx(i)}.cse(2,:),'b');
        ylim([0 1])
        legend('E2','PV2','SOM2','dep','fac')
        title('firing rate(L56)')
        grid on;box on;
    end
    set(gcf,'position',[0 0 800 800])
end

%%=========================================================================
function plot_rates_column2(R)

    figure;
    idx=[2,3,1,4,5];
    for i=1:5

        subplot(3,5,i)% rate(L23)
        hold on;
        plot(R{idx(i)}.gE(4,:),'k','linewidth',2);
        plot(R{idx(i)}.gP(3,:),'c','linewidth',2);
        plot(R{idx(i)}.gS(3,:),'r','linewidth',2);
        plot(R{idx(i)}.cee(4,:),'g');
        plot(R{idx(i)}.cse(4,:),'b');
        ylim([0 1])
        legend('E1','PV1','SOM1','dep','fac')
        title('firing rate(L23 of Col2)')
        grid on;box on;

        subplot(3,5,i+5)% rate(L4)
        hold on;
        plot(R{idx(i)}.gE(6,:),'k','linewidth',2);
        plot(R{idx(i)}.cee(6,:),'g');
        plot(R{idx(i)}.cse(6,:),'b');
        plot(R{idx(i)}.input(1,:),'k');
        ylim([0 1])
        legend('E3','dep','fac','input')
        title('firing rate(L4 of Col2)')
        grid on;box on;

        subplot(3,5,i+10)% rate(L5)
        hold on;
        plot(R{idx(i)}.gE(5,:),'k','linewidth',2);
        plot(R{idx(i)}.gP(4,:),'c','linewidth',2);
        plot(R{idx(i)}.gS(4,:),'r','linewidth',2);
        plot(R{idx(i)}.cee(5,:),'g');
        plot(R{idx(i)}.cse(5,:),'b');
        ylim([0 1])
        legend('E2','PV2','SOM2','dep','fac')
        title('firing rate(L56 of Col2)')
        grid on;box on;
    end
    set(gcf,'position',[0 0 800 800])
end

%==========================================================================
function plot_mua_csd(ref,simMUA,simCSD)
    %--------------
    timeindex=[201:600,1:200,601:1000];
    figure;
    subplot(121);hold on
    h1=plot(ref.MUA(timeindex,:)+repmat(16:-1:1,1000,1),'r','linewidth',2);
    h2=plot(simMUA(timeindex,:)+repmat(16:-1:1,1000,1),'k','linewidth',2);
    set(gca,'ytick',1:16,'yticklabel',16:-1:1,'xtick',200,'xticklabel',200);
    ylim([0 17]);
    plot([1;1]*[200,400,600,800],ylim,'k');
    plot(xlim,[1;1]*(1:16),'k')
    ylabel('channel');xlabel('time (ms)')
    R2=1-sumsqr(ref.MUA(:)-simMUA(:))/sumsqr(ref.MUA(:)-mean(ref.MUA(:))); 
    title(sprintf('MUA (R^2=%g)',R2));
    legend([h1(1),h2(1)],'ref','sim')
    box on;
    %--------------
    subplot(122);hold on
    h1=plot(ref.CSD(timeindex,:)/2+repmat(14:-1:3,1000,1),'r','linewidth',2);
    h2=plot(simCSD(timeindex,:)/2+repmat(14:-1:3,1000,1),'k','linewidth',2);
    set(gca,'ytick',3:14,'yticklabel',14:-1:3,'xtick',200,'xticklabel',200);
    ylim([2 15]);
    plot([1;1]*[200,400,600,800],ylim,'k');
    plot(xlim,[1;1]*(3:14),'k')
    ylabel('channel');xlabel('time (ms)')
    R2=1-sumsqr(ref.CSD(:)-simCSD(:))/sumsqr(ref.CSD(:)-mean(ref.CSD(:))); 
    title(sprintf('CSD (R^2=%g)',R2));
    legend([h1(1),h2(1)],'ref','sim')
    box on;
    set(gcf,'position',[0 0 800 800])

end

%==========================================================================
function plot_mua_csd_detail(R)
    t=1:200;
    figure;
    xtext={'E1','PV1','SOM1','E2','PV2','SOM2','Th','E3'};
    subpl=[2,5,7,3,6,8,1,4];

    xtext2={'E1','E2','E3','PV1','PV2','SOM1','SOM2'};
    nColumns=7;
    for i=1:8

        if i>1 
            subplot(8,nColumns,(i-1)*nColumns+1)
            tmp=[R.E(1:3,:);...
                 R.P(1:2,:);...
                 R.S(1:2,:)];
            plot(t,tmp(i-1,:),'k','linewidth',2);
            ylim([min(tmp(:)) max(tmp(:))]);
            hold on;
            xlim([0 200]);
            plot([1;1]*(51:50:200),ylim,'k');  
            plot(xlim,[0;0],'k'); 
            title(sprintf('PSP (%s)',xtext2{i-1}))  

        end

        subplot(8,nColumns,(subpl(i)-1)*nColumns+2)
            tmp=[R.gE(1,:);...
                 R.gP(1,:);...
                 R.gS(1,:);...
                 R.gE(2,:);...
                 R.gP(2,:);...
                 R.gS(2,:);...
                 R.input;...
                 R.gE(3,:)];
        plot(t,tmp(i,:),'k','linewidth',2);
        xlim([0 200]);ylim([0 ,1]);
        hold on;
        plot([1;1]*(51:50:200),ylim,'k');  
        title(sprintf('firing rate (%s)',xtext{i}))  

        if i>1            
            subplot(8,nColumns,(i-1)*nColumns+3)   
            bar(R.MUAsp(i-1,:));
            xlim([0.5 16.5]);ylim([0,1]);
            set(gca,'xtick',[1,16])
            view([90 90])
            title(sprintf('nMUAsp (%s)',xtext2{i-1})) 
        end
        
        if i>1 
             tmp=[R.input;...
             R.gE(1:3,:);...
             R.gP(1:2,:);...
             R.gS(1:2,:)];
            subplot(8,nColumns,(i-1)*nColumns+4)
            tmp2=R.MUAsp(i-1,:)'*tmp(i,:);
            imagesc(t,1:16,tmp2);
            set(gca,'ytick',[1,16])
            clim([-1 1]*max(abs(tmp2(:))));colormap(flipud(coolwarm));colorbar;
            hold on;plot([1;1]*(51:50:200),ylim,'k');
            title(sprintf('MUA (%s)',xtext2{i-1}))     
        end

        subplot(8,nColumns,(subpl(i)-1)*nColumns+5)
        plot(t,R.flow(:,i),'k','linewidth',2);hold on;
        ylim([0 ,1]*max(R.flow(:)));
        xlim([0 200]);
        plot([1;1]*(51:50:200),ylim,'k'); 
        plot(xlim,[0;0],'k');
        title(sprintf('currents (%s->E)',xtext{i})) 

        subplot(8,nColumns,(subpl(i)-1)*nColumns+6)          
        bar(3:14,R.CSDsp(:,i)/max(abs(R.CSDsp(:))));
        xlim([2.5 14.5]);ylim([-1,1]);
        set(gca,'xtick',[3,14])
        view([90 90])
        title(sprintf('nCSDsp (%s->E)',xtext{i})) 

        subplot(8,nColumns,(subpl(i)-1)*nColumns+7)
        tmp=R.CSDsp(:,i)*R.flow(:,i)';
        imagesc(t,3:14,tmp);
        set(gca,'ytick',[3,14])
        clim([-1 1]*max(abs(tmp(:))));colormap(flipud(coolwarm));colorbar;
        hold on;plot([1;1]*(51:50:200),ylim,'k');
        title(sprintf('CSD (%s->E)',xtext{i}))  
    end
    set(gcf,'position',[0 0 1800 900]);
end
%==========================================================================
function g = model_setting(p,ref)
    %-----inputs -----
    g.nTone       = 1;
    g.soundLength = 2000; % 200 ms
    g.onset       = 0;   
    g.duration    = 2000; % 200 ms
    g.delay       = ref.delayIn*10;
    g.tauin       = 100;  %10 ms decay tau
    %g.ratio      = 0.1;  %decay ratio     
    g.B=0;                % background input   

    %-----network-----                            
    g.Ne=6;   % number of E populations   (E1,E2,E3)
    g.Ns=4;   % number of SOM populations (SOM1,SOM2)
    g.Np=4;   % number of PV populations  (PV1,PV2)
    g.Nt=1;   % number of external input    
    g.w_scaler = 3/10;  % overall scaler
    scaler_e=1;         % from E
    scaler_p=1;         % from PV
    scaler_s=1;         % from SOM
    scaler_L=50;        % lateral inhibition

    %-----
    Pee=[1, 1, 1;...        % {E1,E2,E3}->E1
         1, 1, 1;...        % {E1,E2,E3}->E2  
         1, 1, 1]*p(1);     % {E1,E2,E3}->E3
               
    Ppe=[1, 1, 1;...        % {E1,E2,E3}->PV1
         1, 1, 1]*p(2);     % {E1,E2,E3}->PV2
       
    Pse=[1, 1, 1;...        % {E1,E2,E3}->SOM1
         1, 1, 1]*p(3);     % {E1,E2,E3}->SOM2    

    Pep=[1, 1;...           % {PV1,PV2}->E1 
         1, 1;...           % {PV1,PV2}->E2 
         1, 1]*p(4);        % {PV1,PV2}->E3 
     
    Ppp=[1, 1;...           % {PV1,PV2}->PV1
         1, 1]*p(5);        % {PV1,PV2}->PV2 
     
    Psp=[1, 1;...           % {PV1,PV2}->SOM1
         1, 1]*p(6);        % {PV1,PV2}->SOM2
     
    Pes=[1,1;...            % {SOM1,SOM2}->E1
         1,1;...            % {SOM1,SOM2}->E2
         1,1]*p(7);         % {SOM1,SOM2}->E3
     
    Pps=[1,1;...            % {SOM1,SOM2}->PV1
         1,1]*p(8);         % {SOM1,SOM2}->PV2 
     
    Pss=[0,0;...            % {SOM1,SOM2}->SOM1
         0,0];              % {SOM1,SOM2}->SOM2  

    Pet=[1;...              % extern to E1
         1;...              % extern to E2
         1]*p(9);           % extern to E3
     
    Ppt=[1;...              % extern to PV1
         1]*p(10);          % extern to PV2 
     
    Pst=[0;...              % extern to SOM1
         0];                % extern to SOM2
 
    
    % [ref] Allen brain map
    g.Wee=Pee.*[0.0576, 0.0025, 0.1092;...                      % {E1,E2,E3}->E1
                0.0154, 0.0291, 0.0541;...                      % {E1,E2,E3}->E2  
                0.0054, 0.0007, 0.2017]*g.w_scaler*scaler_e;    % {E1,E2,E3}->E3            

    g.Wpe=Ppe.*[0.3442, 0.0156, 0.3551;...                      % {E1,E2,E3}->PV1
                0.0267, 0.1675, 0.0316]*g.w_scaler*scaler_e;    % {E1,E2,E3}->PV2

    g.Wse=Pse.*[0.1027, 0.0065, 0.2013;...                      % {E1,E2,E3}->SOM1
                0.0135, 0.0264, 0.0166]*g.w_scaler*scaler_e;    % {E1,E2,E3}->SOM2            

    g.Wep=Pep.*[0.1719, 0.0203;...                      % {PV1,PV2}->E1
                0.0092, 0.1387;...                      % {PV1,PV2}->E2
                0.1461, 0.0203]*g.w_scaler*scaler_p;    % {PV1,PV2}->E3

    g.Wpp=Ppp.*[0.1703, 0.0123;...                      % {PV1,PV2}->PV1
                0.0177, 0.1431]*g.w_scaler*scaler_p;    % {PV1,PV2}->PV2           

    g.Wsp=Psp.*[0.0168, 0.0008;...                      % {PV1,PV2}->SOM1
                0.0008, 0.0174]*g.w_scaler*scaler_p;    % {PV1,PV2}->SOM2   

    g.Wes=Pes.*[0.1028, 0.0106;...                      % {SOM1,SOM2}->E1
                0.0015, 0.0322;...                      % {SOM1,SOM2}->E2
                0.0591, 0.0039]*g.w_scaler*scaler_s;    % {SOM1,SOM2}->E3

    g.Wps=Pps.*[0.2268, 0.0008;...                     % {SOM1,SOM2}->PV1
                0.0008, 0.0947]*g.w_scaler*scaler_s;   % {SOM1,SOM2}->PV2

    g.Wss=Pss.*[0.0099, 0.0010;...                     % {SOM1,SOM2}->SOM1
                    0 , 0.0130]*g.w_scaler*scaler_s;   % {SOM1,SOM2}->SOM2

    %[ref] Thalamocortical innervation pattern in mouse auditory and visual cortex (2016)    
    g.Wet=Pet.*[0.225;...            % extern to E1
                0.34 ;...            % extern to E2
                1    ]*g.w_scaler;   % extern to E3
    g.Wpt=Ppt.*[1.25 ;...            % extern to PV1
                1.02 ]*g.w_scaler;   % extern to PV2                     
    g.Wst=Pst.*[0    ;...            % extern to SOM1
                0    ]*g.w_scaler;   % extern to SOM2
            
    %-----two columns-----       
    g.Wee = kron(eye(2),g.Wee);
    g.Wpe = kron(eye(2),g.Wpe);
    g.Wse = kron(eye(2),g.Wse);
    g.Wep = kron(eye(2),g.Wep);
    g.Wpp = kron(eye(2),g.Wpp);
    g.Wsp = kron(eye(2),g.Wsp);
    g.Wes = kron(eye(2),g.Wes);
    g.Wps = kron(eye(2),g.Wps);
    g.Wss = kron(eye(2),g.Wss);
    %----- BF input setting-----
    g.Wse_bf = g.Wse;
    g.Wse_bf(1:2,5) = [0.0065;0.0264]*scaler_L*p(20)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}
    g.Wse_bf(3:4,2) = [0.0065;0.0264]*scaler_L*p(20)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}  
    g.Wet_bf = [g.Wet;g.Wet]*p(29);% BF input
    g.Wpt_bf = [g.Wpt;g.Wpt]*p(29);
    g.Wst_bf = [g.Wst;g.Wst]*p(29);
    %----- nonBF input 1 setting-----
    g.Wse_nbf1 = g.Wse;
    g.Wse_nbf1(1:2,5) = [0.0065;0.0264]*scaler_L*p(21)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}
    g.Wse_nbf1(3:4,2) = [0.0065;0.0264]*scaler_L*p(21)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}
    g.Wet_nbf1 = [g.Wet*p(25);g.Wet]*p(29);% nonBF input
    g.Wpt_nbf1 = [g.Wpt*p(25);g.Wpt]*p(29);
    g.Wst_nbf1 = [g.Wst*p(25);g.Wst]*p(29);
    %----- nonBF input 2 setting-----
    g.Wse_nbf2 = g.Wse;
    g.Wse_nbf2(1:2,5) = [0.0065;0.0264]*scaler_L*p(22)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}
    g.Wse_nbf2(3:4,2) = [0.0065;0.0264]*scaler_L*p(22)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}
    g.Wet_nbf2 = [g.Wet*p(26);g.Wet]*p(29);% nonBF input
    g.Wpt_nbf2 = [g.Wpt*p(26);g.Wpt]*p(29);
    g.Wst_nbf2 = [g.Wst*p(26);g.Wst]*p(29);
    %----- nonBF input 3 setting-----
    g.Wse_nbf3 = g.Wse;
    g.Wse_nbf3(1:2,5) = [0.0065;0.0264]*scaler_L*p(23)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}
    g.Wse_nbf3(3:4,2) = [0.0065;0.0264]*scaler_L*p(23)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2} 
    g.Wet_nbf3 = [g.Wet*p(27);g.Wet]*p(29);% nonBF input
    g.Wpt_nbf3 = [g.Wpt*p(27);g.Wpt]*p(29);
    g.Wst_nbf3 = [g.Wst*p(27);g.Wst]*p(29);
    %----- nonBF input 4 setting-----
    g.Wse_nbf4 = g.Wse;
    g.Wse_nbf4(1:2,5) = [0.0065;0.0264]*scaler_L*p(24)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}
    g.Wse_nbf4(3:4,2) = [0.0065;0.0264]*scaler_L*p(24)*g.w_scaler*scaler_e; % lateral inhibition: E2  -> {S1,S2}
    g.Wet_nbf4 = [g.Wet*p(28);g.Wet]*p(29);% nonBF input
    g.Wpt_nbf4 = [g.Wpt*p(28);g.Wpt]*p(29);
    g.Wst_nbf4 = [g.Wst*p(28);g.Wst]*p(29);
        
    %-----short-term plasticity -----  
    g.STPee=1;                % depression on Wee
    g.tauree=200*1e-3;        % time constant
    g.kree=20*p(11);          % depressing rate 
    
    g.STPse=1;                % facilitation on Wse
    g.Use=0.05;               % initial U 
    g.taufse=670*1e-3;        % time constant
    g.kfse=600*p(12);         % facilitating rate 
    
    g.STPth=0;                % depression on Wet,Wst,Wpt (from thalamus)
    g.taurth=90*1e-3;         % time constant
    g.krth=5*59.4;            % depressing rate 
    
    %-----synaptic kernel time constant-----
    % [ref] Data-driven multiscale model of macaque auditory thalamocortical 
    %       circuits reproduces in vivo dynamics, bioRxiv, 2022 
    %       E->E:   AMPA       (rise/fall = 1/5.3 ms) 
    %               NMDA       (rise/fall = 3/70 ms) 
    %       SOM->E: slow GABAA (rise/fall = 2/100 ms)
    %               GABAB      (rise/fall = 25/300 ms)
    scaler_tau = p(13);   

    %---rise tau---
    g.taue1=scaler_tau*[1,4.5,2.1]*1e-3; % from E to (E(AMPA),S,P) 
    g.taus1=scaler_tau*[2,2]*1e-3;       % from SOM to (E,P) 
    g.taup1=scaler_tau*[1,1.4,3.5]*1e-3; % from PV to (E,S,P)
    g.taueNMDA1=scaler_tau*[3]*1e-3;     % E-E (NMDA)
    g.tausGABAB1=scaler_tau*[25]*1e-3;   % SOM->E(GABA_B)
    %---fall tau---
    g.taue2=scaler_tau*[5.3,25.2,5.6]*1e-3; % from E to (E,S,P)
    g.taus2=scaler_tau*[100,100]*1e-3;      % from SOM to (E,P) 
    g.taup2=scaler_tau*[18.2,101,5.5]*1e-3; % from PV to (E,S,P)
    g.taueNMDA2=scaler_tau*[70]*1e-3;       % E-E (NMDA)
    g.tausGABAB2=scaler_tau*[300]*1e-3;     % SOM->E(GABA_B)

    g.He=[14400,3090,7250]/scaler_tau;      % from E to (E,S,P) 
    g.Hs0=[-1800,-1800]/scaler_tau;         % from SOM to (E,P)  
    g.Hp0=[-4000,-7380,-5530]/scaler_tau;   % from PV to (E,S,P)  
    g.HeNMDA=1200/scaler_tau;
    g.HsGABAB0=-100/scaler_tau;
    
    g.Hs=abs(g.Hs0); g.Hp=abs(g.Hp0);
    g.HsGABAB=abs(g.HsGABAB0);           
   
    g.AMPA_ratio=0.83; % E->E connection
    g.GABAB_ratio=0.5; % E->SOM connection
    g.Th_NMDA_ratio=1; % 0~1
    %-----sigmoid function------
    scaler_slope = p(14);
    g.theta_e            = scaler_slope*6   *1e-3;        % E firing 50%  
    g.theta_s            = scaler_slope*2.76*1e-3;        % SOM firing 50% 
    g.theta_p            = scaler_slope*15.6*1e-3;        % PV firing 50% 
    g.slope_e            = 620/scaler_slope;              % E slope (pps/V)   
    g.slope_s            = 1140/scaler_slope;             % SOM slope (pps/V) 
    g.slope_p            = 290/scaler_slope;              % PV slope (pps/V)    
end
   
   

