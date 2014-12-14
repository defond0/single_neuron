close all
clear all
% An implementation of Eugene M. Izhikevich's Simple Model of Spiking
% Neurons 2003 fit into a Linear Hybrid Automota as described in Symbolic Analysis of the Neuron
% Action Potential by Pei Ye∗, Emilia Entcheva†, Scott A. Smolka∗, Radu Grosu

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

iv = c;
iu = b*v;

resting = 'resting';
stimulated = 'stimulated';
rising = 'rising';
falling = 'falling';

time = 1000;
spikeonset=250
spikeduration=5;

spikeonset1=500
spikeduration1=5;

voltages=zeros(time,1);
Ivoltages = zeros(time,1);
Mode = resting

I=0; %intial input of zero
for t=1:time %simulate for 1 second
   if t==spikeonset % send dc pulse of I=10 
       I=10;
   end
   
   if t==spikeonset+spikeduration;
       I=0;
   end
   iv = iv + ((.04*(iv^2))+(5*iv) +140 -iu + I);
   iu = iu + a*((b*iv) - iu ); 
   Ivoltages(t)= iv;
   if(iv>=30)
       iv=c;
       iu = iu+d;
   end
      mini = min(Ivoltages);
      maxi = max(Ivoltages);
   
    
   q=.01;
  
   adj=77;
      
   switch Mode
       case resting
           if(I==0)
                v = v  -(.9*v) -u -adj;
                u = u + a*((b*v) - u );
           else
               v = v + abs((q*v) -u -adj + I);
               u = u + a*((b*v) - u );
               Mode=stimulated
           end
       case stimulated
               v = v + abs((q*v) -u -adj + I);
               u = u + a*((b*v) - u );
               if((I==0)&(v<30))
                   Mode=resting
               elseif(v>=30)
                   Mode=rising
               end
       case rising
                if (v>=90)
                     v=c;
                     u = u+d;
                     Mode=falling
                end
                v =v+ abs((q*v) -u -adj);
                u = u + a*((b*v) - u );
                
       case falling
            v = (v  -(.9*v) -u)-adj;
            u = u + a*((b*v) - u );
            if(v>=c-10)
                Mode=resting
                
            end
   end

   voltages(t)= v;
   
end

dif = abs(voltages-Ivoltages);
averageError = sum(dif)/time
figure
plot(voltages,'b');
hold on;
plot(Ivoltages,'r');
hold on;
plot(dif,'g');
hold off;
xlabel('Time Ms')
ylabel('Membrane Potential mV')
axis([0 time -100 100])


