
% An implementation of Eugene M. Izhikevich's Simple Model of Spiking Neurons 2003
close all
clear all


neuron_type = 'rs' %Excititory (rs, ib, ch) Inhibitory (fs, lts, tc, rz) 


switch neuron_type
    case 'rs'
        a = .02; %scale of recovery variable, smaller value, slower recovery
        b = .2; %sensitivity of the recovery variable u to subthreshold fluctuations
        c = -65; %describes resting potential
        d = 8; %afterspike jump
    case 'ib'
        a = .02; %scale of recovery variable, smaller value, slower recovery
        b = .2; %sensitivity of the recovery variable u to subthreshold fluctuations
        c = -55; %describes resting potential
        d = 4; %afterspike jump
    case 'ch'
        a = .02; %scale of recovery variable, smaller value, slower recovery
        b = .2; %sensitivity of the recovery variable u to subthreshold fluctuations
        c = -50; %describes resting potential
        d = 2; %afterspike jump
    case 'fs'
        a = .1; %scale of recovery variable, smaller value, slower recovery
        b = .2; %sensitivity of the recovery variable u to subthreshold fluctuations
        c = -65; %describes resting potential
        d = 2; %afterspike jump
    case 'lts'
        a = .02; %scale of recovery variable, smaller value, slower recovery
        b = .25; %sensitivity of the recovery variable u to subthreshold fluctuations
        c = -65; %describes resting potential
        d = 2; %afterspike jump
    case 'tc'
        a = .02; %scale of recovery variable, smaller value, slower recovery
        b = .2; %sensitivity of the recovery variable u to subthreshold fluctuations
        c = -60; %describes resting potential
        d = 2; %afterspike jump
    case 'rz'
        a = .1; %scale of recovery variable, smaller value, slower recovery
        b = .26; %sensitivity of the recovery variable u to subthreshold fluctuations
        c = -65; %describes resting potential
        d = 2; %afterspike jump
end

v = c ;
u = b*v;

voltages=zeros(1000,1);

I=0; %intial input of zero
for t=1:1000 %simulate for 1 second
   if(v>=30) % reset after spike
       v = c; % reset voltage to rest
       u = u + d; %reset recovery variable
   end
   if t==100 % send dc pulse of I=10 at ten seconds, see paper
       I=10;
   end
   if strcmp(neuron_type,'rz') && t==200
       I=I+10;
   end 
   if strcmp(neuron_type,'rz') && t==250
       I=I-10;
   end 
   v = v + ((.04*(v^2))+(5*v) +140 -u + I);
   u = u + a*((b*v) - u ); 
   voltages(t)= v;
   
end

figure
plot(voltages)
xlabel('Time Ms')
ylabel('Membrane Potential mV')
axis([0 1000 -100 30])



