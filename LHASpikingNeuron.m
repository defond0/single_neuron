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
        restingOut=zeros(3,1);
        stimulatedOut=zeros(3,1);
        risingOut=zeros(3,1);
        fallingOut=zeros(3,1);
        stimThresh=-30;
        fallingThresh=-8;%25
        risingThresh=-10;
    end
    
    methods
        %output function to get the lha description after analysis
        %has been run
        function [restingOut,stimulatedOut,risingOut,fallingOut]=getLHA(LHASpikingNeuron);
            restingOut=LHASpikingNeuron.restingOut;
             stimulatedOut=LHASpikingNeuron.stimulatedOut;
             risingOut=LHASpikingNeuron.risingOut
             fallingOut=LHASpikingNeuron.fallingOut
        end
        
        %construct the LHASpikingNeuron object with simulation time, onset
        %of stimulation, duration of stimulation. v is the 2nd order model,
        %vdt is the derivative wit respect to time of this model, and fb is the feedback relation. 
        function lha = LHASpikingNeuron(time,onset,duration,v,vdt,fb)
            lha.time=time;
            lha.onset=onset;
            lha.duration=duration;
            lha.v=v;
            lha.vdt=vdt;
            lha.fb=fb;
        end     
        
        %Try values of A that have the lowest absolute error. Mode is which
        %mode is to be solved for, index refers to what indexes on the
        %given model's simulation the model was in mode. voltages are the
        %voltages from running the given model.
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
        
        %Analogous to VA, accept now we find a constant term for the given
        %mode and agiven A. 
        function Cee = findC(LHASpikingNeuron,A,mode,index,voltages,Avoltages)
            time = LHASpikingNeuron.time;
            onset = LHASpikingNeuron.onset;
            duration = LHASpikingNeuron.duration;
            
            currents=zeros(time,1);
            currents([onset:onset+duration])=10;
            c = -65;
            d = 8;
            C = 0;
            prevC =0;
            p=LHASpikingNeuron.fb(1);
            q=LHASpikingNeuron.fb(2);
          
            stimulated = 'stimulated';
            
            
            stimThresh=LHASpikingNeuron.stimThresh;
            fallingThresh=LHASpikingNeuron.fallingThresh;
            risingThresh=LHASpikingNeuron.risingThresh;
            
            found=false;
            prevAverageDif=inf; 
            while(~found)  
            v=c;
            u= v*p/-q;        
            Ivoltages = zeros(time,1);
            
          
            for t=1:time 
                    I=currents(t);
                     testStim=~(sum(find(index==t))==0);
                     if(testStim)
                        v=v+v*A+C;
                        if(strcmp(mode,stimulated))
                            v=v+I;
                        end
                     else
                        v = Avoltages(t); 
                     end
                    Ivoltages(t)=v;
            end
            dif = abs(voltages-Ivoltages);
            averageDiff = sum(dif)/time;
            
            if(prevAverageDif<=averageDiff)
                    averageDiff=prevAverageDif                    
                    C=prevC;
                    found=true;
                 else
                    prevAverageDif=averageDiff;                   
                    prevC=C;
                    C=C+1;
                 end  
            end
            Cee=C
        end
        
        
        %A method to generate Avoltages for findC
        function [IV,D]= prepForC(LHASpikingNeuron,restA,stimA,risingA,fallingA,voltages) 
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
                       if(v>=stimThresh)
                           Mode=rising
                       
                       elseif(I==0)&(v<stimThresh)
                           v = c-d;
                           u = u + d;
                           Mode=resting
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
            IV = Ivoltages;
            D = dif;
           
          
        end
        
        %A method to generate the final lha and the difference between it
        %and the model
        function [IV,D]= finalizeAndCheck(LHASpikingNeuron,restA,restC,stimA,stimC,risingA,risingC,fallingA,fallingC,voltages) 
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
                        v=v+(v*stimA)+I+stimC;
                        Mode=stimulated
                      else
                          v = v + v*restA+restC;   
                      end                         
                   case stimulated        
                       v = v + v*stimA+I+stimC;
                       if(v>=stimThresh)
                           Mode=rising
                       
                       elseif(I==0)&(v<stimThresh)
                           v = c-d;
                           u = u + d;
                           Mode=resting
                       end

                   case rising
                        if(v>risingThresh)
                            v = c-20;
                            u = u + d;
                            Mode=falling
                        else
                            v = v + v*risingA+I+risingC;
                       end

                   case falling
                       v=v+fallingA*v+fallingC;   
                       if(v>=c+fallingThresh)
                           Mode=resting
                  
                       end
                end
                Ivoltages(t)=v;
            end           
            dif = abs(voltages-Ivoltages);
            averageDiff = sum(dif)/time
            IV = Ivoltages;
            D = dif;
           
          
        end
        
        %get voltages for the given model
        function voltagesOut = runGivenModel(LHASpikingNeuron)
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
        
        %generate the output for the lha. This method makes calls to ones
        %above to get an equation of the form V= V+A*V'+C for all modes,
        %risng, stimulated, falling and rising. 
        function [IV,V,restingOut,stimulatedOut,risingOut,fallingOut]=lha(LHASpikingNeuron)
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
            stimA=LHASpikingNeuron.findVA(stimulated,stimulatedIdx,voltages)
            risingA=LHASpikingNeuron.findVA(rising,risingIdx,voltages)
            fallingA=LHASpikingNeuron.findVA(falling,fallingIdx,voltages)

          
            [Ivoltages,dif]=LHASpikingNeuron.prepForC(restA,stimA,risingA,fallingA,voltages) ;
          
            
            restC=LHASpikingNeuron.findC(restA,resting,restingIdx,voltages,Ivoltages);
            stimC=LHASpikingNeuron.findC(stimA,stimulated,stimulatedIdx,voltages,Ivoltages);
            risingC=LHASpikingNeuron.findC(risingA,rising,risingIdx,voltages,Ivoltages);
            fallingC=LHASpikingNeuron.findC(fallingA,falling,fallingIdx,voltages,Ivoltages);
            
            [Ivoltages,dif]=LHASpikingNeuron.finalizeAndCheck(restA,restC,stimA,stimC,risingA,risingC,fallingA,fallingC,voltages) ;
            
            restingOut=zeros(3,1);
            stimulatedOut=zeros(3,1);
            risingOut=zeros(3,1);
            fallingOut=zeros(3,1);
            
            restingOut(2)=restA;
            stimulatedOut(2)=stimA;
            risingOut(2)=risingA;
            fallingOut(2)=fallingA;
            
            
            restingOut(3)=restC;
            stimulatedOut(3)=stimC;
            risingOut(3)=risingC;
            fallingOut(3)=fallingC;
            
            
            LHASpiking.restingOut=restingOut        ;
            LHASpiking.stimulatedOut=stimulatedOut;
            LHASpiking.risingOut=risingOut;
            LHASpiking.fallingOut=fallingOut;
            
            IV=Ivoltages;
            V=voltages;
         
            
            
        end
        
        
        
    end
end

