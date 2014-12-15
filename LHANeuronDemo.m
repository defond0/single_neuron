

neuron_type = 'fs' %Excititory (rs, ib, ch) Inhibitory (fs, lts, tc, rz) 

duration = 100;

switch neuron_type
    case 'rs'
        a = .02; %scale of recovery variable, smaller value, slower recovery
        b = .2; %sensitivity of the recovery variable u to subthreshold fluctuations
        c = -65; %describes resting potential
        d = 8; %afterspike jump
    case 'ib'
        a = .02;
        b = .2; 
        c = -55; 
        d = 4; 
    case 'ch'
        a = .02; 
        b = .2; 
        c = -50; 
        d = 2; 
    case 'fs'
        a = .1; 
        b = .2; 
        c = -65; 
        d = 2; 
    case 'lts'
        a = .02;
        b = .25; 
        c = -65; 
        d = 2; 
    case 'tc'
        a = .02; 
        b = .2; 
        c = -60;
        d = 2; 
    case 'rz'
        a = .1;
        b = .26; 
        c = -65; 
        d = 2; 
end

%2nd order, 1st order, constant
v =[.04,5,140];
vdt = [0,.08,0];
%v term, u term
fb = [a*b,-a];

lha = LHASpikingNeuron(1000,250,duration,v,vdt,fb);
lha.lha();




     
 