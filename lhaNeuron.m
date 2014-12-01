close all
clear all

% An implementation of Eugene M. Izhikevich's Simple Model of Spiking
% Neurons 2003 fit into a Linear Hybrid Automota as described in Symbolic Analysis of the Neuron
% Action Potential by Pei Ye∗, Emilia Entcheva†, Scott A. Smolka∗, Radu Grosu

neuron_type = 'ib' %Excititory (rs, ib, ch) Inhibitory (fs, lts, tc, rz) 


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

f = 'falling';
s = 'stimulated';
r = 'rising';

voltages=zeros(2000,1);

Mode = f

I=0; %intial input of zero
for t=1:2000 %simulate for 1 second
    
   if t==100 % send dc pulse of I=10 at ten seconds, see paper
       I=10;
   end
   
   if t==1000
       I=0;
   end
    
   switch Mode
       case f
           v = v + ((.04*(v^2))+(5*v) +140 -u);
           u = u + a*((b*v) - u );
           if(I>0)
               Mode = s
           end
       case s
           v = v + ((.04*(v^2))+(5*v) +140 -u + I);
           u = u + a*((b*v) - u );
           if(v>=30) % reset after spike
                v = c; % reset voltage to rest
                u = u + d; %reset recovery variable
                Mode = r
           end
           if(I==0)
               Mode = f 
           end
       case 'rising'
           v = v + ((.04*(v^2))+(5*v) +140 -u + I);
           u = u + a*((b*v) - u ); 
           if(v>=30) % reset after spike
                v = c; % reset voltage to rest
                u = u + d; %reset recovery variable
                Mode = f
           end
   end
   voltages(t)= v;
   
end

figure
plot(voltages)
xlabel('Time Ms')
ylabel('Membrane Potential mV')
axis([0 2000 -100 30])

