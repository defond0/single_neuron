
% An implementation of Eugene M. Izhikevich's Simple Model of Spiking Neurons 2003



neuron_type = 'rs' %Excititory (rs, ib, ch) Inhibitory (fs, lts, tc, rz) 


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

v = c ;
u = b*v;

time =1000

voltages=zeros(time,1);

I=0; %intial input of zero
for t=1:time %simulate for 2 second
   if t==250 % send dc pulse of I=10 at ten seconds, see paper
       I=10;
   end
   
   if t==500
       I=0;
   end
   if(v>=30) % reset after spike
       v = c; % reset voltage to rest
       u = u + d; %reset recovery variable
   end
   v = v + ((.04*(v^2))+(5*v) +140 -u + I);
   u = u + a*((b*v) - u ); 
   voltages(t)= v;
   
end

figure
plot(voltages)
xlabel('Time Ms')
ylabel('Membrane Potential mV')
axis([0 time -100 30])



