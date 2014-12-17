classdef LHASpikingNeuron
    %Generate Linear Hybrid Automota from 2nd order neuron model
    %   
    
    properties
        time;
        onset;
        duration;
        v;
        vdt;
        fb; 
        restingOut;
        stimulatedOut;
        risingOut;
        fallingOut;
        stimThresh=0;
        fallingThresh=-8;%25
        risingThresh=-80;
    end
    
    methods
        
        function lha = LHASpikingNeuron(time,onset,duration,v,vdt,fb)
            lha.time=time;
            lha.onset=onset;
            lha.duration=duration;
            lha.v=v;
            lha.vdt=vdt;
            lha.fb=fb;
        end     
        
        function VA = findVA(LHASpikingNeuron,mode,index,voltages)
            c = -65;
            d = 8;
            p=LHASpikingNeuron.fb(1);
            q=LHASpikingNeuron.fb(2);
            
           
            time = LHASpikingNeuron.time;
            onset = LHASpikingNeuron.onset;
            duration = LHASpikingNeuron.duration;
            
            i=0;
            prevAverageDif=inf;
            found=false;
            A=LHASpikingNeuron.vdt(2);
            prevA=0;
            stimulated = 'stimulated';
            currents=zeros(time,1);
            currents([onset:onset+duration])=10;
            while(~found)            
               Ivoltages = zeros(time,1);
                v=c;                
                u= v*p/-q;
                for t=1:time 
                    I=currents(t);
                     testStim=~(sum(find(index==t))==0);
                     if(testStim)
                        v=v+v*A;
                        if(strcmp(mode,stimulated))
                            v=v+I;
                        end
                     else
                        v = voltages(t); 
                     end
                    Ivoltages(t)=v;
                end               
                dif = abs(voltages-Ivoltages);               
                averageDiff = sum(dif)/time;
                if(prevAverageDif<=averageDiff)
                    averageDiff=prevAverageDif                    
                    A=prevA;
                    found=true;
                 else
                    prevAverageDif=averageDiff;                   
                    prevA=A;
                    A=abs(A)*-.1;
                 end  
                 i=i+1;
            end
            VA=A;
        end
        
        function VA = findAK(LHASpikingNeuron,index,voltages,cluster)
            c = -65;
            d = 8;
            p=LHASpikingNeuron.fb(1);
            q=LHASpikingNeuron.fb(2);
                     
            time = LHASpikingNeuron.time;
            onset = LHASpikingNeuron.onset;
            duration = LHASpikingNeuron.duration;
            
            i=0;
            prevAverageDif=inf;
            found=false;
            A=LHASpikingNeuron.vdt(2);
            prevA=0;
            currents=zeros(time,1);
            currents([onset:onset+duration])=10;
            while(~found)            
               Ivoltages = zeros(time,1);
                v=c;                
                u= v*p/-q;
                for t=1:time 
                    I=currents(t);
                     if(index(t)==cluster)
                        v=v+v*A+I;
                     else
                        v = voltages(t); 
                     end
                    Ivoltages(t)=v;
                end               
                dif = abs(voltages-Ivoltages);               
                averageDiff = sum(dif)/time;
                if(prevAverageDif<=averageDiff)
                    averageDiff=prevAverageDif                    
                    A=prevA;
                    found=true;
                 else
                    prevAverageDif=averageDiff;                   
                    prevA=A;
                    A=abs(A)*-.1;
                 end  
                 i=i+1;
            end
            VA=A;
        end
        
        function finalizeAndCheck(LHASpikingNeuron,restA,stimA,risingA,fallingA, stimThresh, fallingThresh,risingThresh,voltages) 
            time = LHASpikingNeuron.time;
            onset = LHASpikingNeuron.onset;
            duration = LHASpikingNeuron.duration;
            
            currents=zeros(time,1);
            currents([onset:onset+duration])=10;
            
            c = -65;
            d = 8;
            
            resting = 'resting';
            stimulated = 'stimulated';
            rising = 'rising';
            falling = 'falling';
            
            stimThresh=LHASpikingNeuron.stimThresh;
            fallingThresh=LHASpikingNeuron.fallingThresh;
            risingThresh=LHASpikingNeuron.risingThresh;
            Ivoltages = zeros(time,1);
            p=LHASpikingNeuron.fb(1);
            q=LHASpikingNeuron.fb(2);
            v=c;
            u= v*p/-q;        
            Mode=resting
            for t=1:time 
                I=currents(t);                                 
                switch Mode
                   case resting                                                                              
                      if(I~=0)
                        v=v+(v*stimA)+I;
                        Mode=stimulated
                      else
                          v = v + v*restA;   
                      end                         
                   case stimulated        
                       v = v + v*stimA+I;                                                 
                       if(I==0)&(v<stimThresh)
                           v = c-d;
                           u = u + d;
                           Mode=resting
                       elseif(v>stimThresh)
                           Mode=rising
                       end

                   case rising
                        if(v>risingThresh)
                            v = c-20;
                            u = u + d;
                            Mode=falling
                        else
                            v = v + v*risingA+I;
                       end

                   case falling
                       v=v+fallingA*v;   
                       if(v>=c+fallingThresh)
                           Mode=resting
                  
                       end
                end
                Ivoltages(t)=v;
            end           
            dif = abs(voltages-Ivoltages);
            averageDiff = sum(dif)/time
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
          
        end
        
        function voltagesOut = runGivenModel(LHASpikingNeuron);
            vsq=LHASpikingNeuron.v(1);
            vl=LHASpikingNeuron.v(2);
            cons=LHASpikingNeuron.v(3);
     
            p=LHASpikingNeuron.fb(1);
            q=LHASpikingNeuron.fb(2);
            
            time = LHASpikingNeuron.time;
            onset = LHASpikingNeuron.onset;
            duration = LHASpikingNeuron.duration;
            voltages=zeros(time,1);
            currents=zeros(time,1);
            currents([onset:onset+duration])=10;
            
            c = -65;
            d = 8;
            v=c;
            u= v*p/-q;
            for t=1:time 
                I=currents(t);              
                if(v>=30) 
                   v = c;
                   u = u + d; 
                end
                 v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                 u = u + (p*v + q*u);
                 voltages(t)=v;
             end
             voltagesOut=voltages;
            
           
        end
          
        function  [idx,C,sumd,D]=lhaKmean(LHASpikingNeuron)            
            voltages = LHASpikingNeuron.runGivenModel()
            time = LHASpikingNeuron.time;
            [idx,C,sumd,D] = kmeans(voltages,3);
            [m, maxIdx] = max(C);
            [m, minIdx] = max(C);
            %cluster1A = LHASpikingNeuron.findAK(idx,voltages,1)
            %cluster2A = LHASpikingNeuron.findAK(idx,voltages,2)
            %cluster3A = LHASpikingNeuron.findAK(idx,voltages,3)
            
            
            
        end
        
        function voltages=lha(LHASpikingNeuron)
            c = -65;
            d = 8;
            
            stimThresh=LHASpikingNeuron.stimThresh;
            fallingThresh=LHASpikingNeuron.fallingThresh;
            risingThresh=LHASpikingNeuron.risingThresh;
            
            vsq=LHASpikingNeuron.v(1);
            vl=LHASpikingNeuron.v(2);
            cons=LHASpikingNeuron.v(3);
     
            p=LHASpikingNeuron.fb(1);
            q=LHASpikingNeuron.fb(2);
            
            time = LHASpikingNeuron.time;
            onset = LHASpikingNeuron.onset;
            duration = LHASpikingNeuron.duration;
            
            resting = 'resting';
            stimulated = 'stimulated';
            rising = 'rising';
            falling = 'falling';
            
            v=c;
            u= v*p/-q;
            
            Mode=resting
            
            voltages=zeros(time,1);
           
            restV=zeros(time,1);
            stimulatedV=zeros(time,1);
            risingV=zeros(time,1);
            fallingV=zeros(time,1);
            currents=zeros(time,1);
            currents([onset:onset+duration])=10;
            
            for t=1:time 
                I=currents(t);              
                if(v>=30) 
                   v = c;
                   u = u + d; 
                end
                 v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                 u = u + (p*v + q*u); 
                 switch Mode
                   case resting                              
                       if(I~=0)
                           Mode=stimulated
                           stimV(t)=v;
                       else
                            restV(t)=v;
                       end
                   case stimulated
                       stimulatedV(t)=v ; 
                       if(I==0)&(v<stimThresh)
                           Mode=resting
                       elseif(v>stimThresh)
                           Mode=rising
                       end
                   case rising
                       risingV(t)=v;
                        if(v>risingThresh)
                             Mode=falling;
                            
                      
                       end
                   case falling
                       fallingV(t)=v;
                       if(v>=c+fallingThresh)
                           Mode=resting
                       end
                 end                 
                 voltages(t)= v;                  
            end
   
            
            restingIdx=find(restV);
            stimulatedIdx=find(stimulatedV);
            risingIdx=find(risingV);
            fallingIdx=find(fallingV);
            
           
            
            
            restA=LHASpikingNeuron.findVA(resting,restingIdx,voltages)
            %pause = input('finished with restA');
            stimA=LHASpikingNeuron.findVA(stimulated,stimulatedIdx,voltages)
            %pause = input('finished with stimA');
            risingA=LHASpikingNeuron.findVA(rising,risingIdx,voltages)
            %pause = input('finished with risingA');
            fallingA=LHASpikingNeuron.findVA(falling,fallingIdx,voltages)
            %pause = input('finished with fallingA');
            
            LHASpikingNeuron.finalizeAndCheck(restA,stimA,risingA,fallingA,stimThresh, fallingThresh,risingThresh,voltages) ;
            
            restingOut=zeros(3,1);
            stimulatedOut=zeros(3,1);
            risingOut=zeros(3,1);
            fallingOut=zeros(3,1);
            
            restingOut(2)=restA;
            stimulatedOut(2)=stimA;
            risingOut(2)=risingA;
            fallingOut(2)=fallingA;
            
            LHASpiking.restingOut=restingOut   ;         
            LHASpiking.stimulatedOut=stimulatedOut;
            LHASpiking.risingOut=risingOut;
            LHASpiking.fallingOut=fallingOut;
            
            
            
            
        end
        
        
        
    end
end

