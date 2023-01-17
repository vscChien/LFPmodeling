% Output:
%      [to x from]
%  Wee {E1,E2,E3}x{E1,E2,E3}
%  Wpe {PV1,PV2}x{E1,E2,E3}
%  Wse {SOM1,SOM2}x{E1,E2,E3}
%  Wep {E1,E2,E3}x{PV1,PV2}
%  Wpp {PV1,PV2}x{PV1,PV2}
%  Wsp {SOM1,SOM2}x{PV1,PV2}
%  Wes {E1,E2,E3}x{SOM1,SOM2}
%  Wps {PV1,PV2}x{SOM1,SOM2}
%  Wss {SOM1,SOM2}x{SOM1,SOM2}
% ,where 
%  E1: E2/3
%  E2: E5/6
%  E3: E4
%  PV1: PV2/3/4
%  PV2: PV5/6
%  SOM1: SOM2/3/4
%  SOM2: SOM5/6
%
% [ref]
% Billeh, Yazan N., et al. "Systematic integration of structural and 
% functional data into multi-scale models of mouse primary visual cortex." 
% Neuron 106.3 (2020): 388-403.
function [Wee,Wpe,Wse,Wep,Wpp,Wsp,Wes,Wps,Wss]=get_connectivity()
    
    % Connection probabilities at 75um
    % E2/3, PV2/3, SST2/3, E4, PV4, SST4, E5, PV5, SST5, E6, PV6, SST6
    connP=[0.16,0.395,0.182,0.016,0.083,0.083,0.083,0.081,0.102,0,0,0;...
    0.411,0.451,0.03,0.05,0.05,0.05,0.07,0.073,0,0,0,0;...
    0.424,0.857,0.082,0.05,0.05,0.05,0.021,0,0,0,0,0;...
    0.14,0.1,0.1,0.243,0.43,0.571,0.104,0.101,0.128,0.032,0,0;...
    0.25,0.05,0.05,0.437,0.451,0.03,0.088,0.091,0.03,0,0,0;...
    0.25,0.05,0.05,0.351,0.857,0.082,0.026,0.03,0,0,0,0;...
    0.021,0.05,0.05,0.007,0.05,0.05,0.116,0.083,0.063,0.047,0.03,0.03;...
    0,0.102,0,0,0.034,0.03,0.455,0.361,0.03,0.03,0.01,0.01;...
    0.169,0,0.017,0.056,0.03,0.006,0.317,0.857,0.04,0.03,0.01,0.01;...
    0,0,0,0,0,0,0.012,0.01,0.01,0.026,0.145,0.1;...
    0.1,0,0,0.1,0,0,0.1,0.03,0.03,0.1,0.08,0.1;...
    0,0,0,0,0,0,0.03,0.03,0.03,0.1,0.05,0.05]'; % to x from
       
    % Sznaptic strengths, unitary PSP (mV)
    % E2/3, PV2/3, SST2/3, E4, PV4, SST4, E5, PV5, SST5, E6, PV6, SST6
    connS=[0.36,1.49,0.86,0.34,1.39,0.69,0.74,1.32,0.53,0,0,0;...
    0.48,0.68,0.42,0.56,0.68,0.42,0.2,0.79,0,0,0,0;...
    0.31,0.5,0.15,0.3,0.5,0.15,0.22,0,0,0,0,0;...
    0.78,1.39,0.69,0.83,1.29,0.51,0.63,1.25,0.52,0.96,0,0;...
    0.56,0.68,0.42,0.64,0.68,0.42,0.73,0.94,0.42,0,0,0;...
    0.3,0.5,0.15,0.29,0.5,0.15,0.28,0.45,0.28,0,0,0;...
    0.47,1.25,0.52,0.38,1.25,0.52,0.75,1.2,0.52,0.4,2.5,0.52;...
    0,0.51,0,0,0.94,0.42,0.81,1.19,0.41,0.81,1.19,0.41;...
    0.25,0,0.39,0.28,0.45,0.28,0.27,0.4,0.4,0.27,0.4,0.4;...
    0,0,0,0,0,0,0.23,2.5,0.52,0.94,3.8,0.52;...
    0.81,0,0,0.81,0,0,0.81,1.19,0.41,0.81,1.19,0.41;...
    0,0,0,0,0,0,0.27,0.4,0.4,0.27,0.4,0.4]'; % to x from
   
    
    % merged_conn = mean(prob) * mean(strength)
    % E  : L23,L4,L56
    % PV : L234,L56
    % SOM: L234,L56
    connP=[connP(:,1),(connP(:,2)+connP(:,5))/2,(connP(:,3)+connP(:,6))/2,connP(:,[4 7 8 9])]; % mean merge PV,SOM in L23 and L4
    connP=[connP(1,:);(connP(2,:)+connP(5,:))/2;(connP(3,:)+connP(6,:))/2;connP([4 7 8 9],:)]; % mean merge PV,SOM in L23 and L4
    connS=[connS(:,1),(connS(:,2)+connS(:,5))/2,(connS(:,3)+connS(:,6))/2,connS(:,[4 7 8 9])]; % mean merge PV,SOM in L23 and L4
    connS=[connS(1,:);(connS(2,:)+connS(5,:))/2;(connS(3,:)+connS(6,:))/2;connS([4 7 8 9],:)]; % mean merge PV,SOM in L23 and L4
    conn=connP.*connS;
    
    iE=[1,5,4];%E1,E2,E3
    iP=[2,6];  %PV1,PV2
    iS=[3,7];  %SOM1,SOM2
    
    Wee=conn(iE,iE);
    Wpe=conn(iP,iE);
    Wse=conn(iS,iE);
    
    Wep=conn(iE,iP);
    Wpp=conn(iP,iP);
    Wsp=conn(iS,iP);
    
    Wes=conn(iE,iS);
    Wps=conn(iP,iS);
    Wss=conn(iS,iS);
    
end
