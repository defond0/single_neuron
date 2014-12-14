classdef LHANeuron
    %LHANeuron this model can generate a LinearHybrid Automota for a 2nd order model of a
    %neuron. 
    %   Assumptions:
    
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
    end
    
    methods 
        function lha = LHANeuron(time,onset,duration,v,vdt,fb)
            lha.time=time;
            lha.onset=onset;
            lha.duration=duration;
            lha.v=v;
            lha.vdt=vdt;
            lha.fb=fb;
        end
        
        function averageError = lha(LHANeuron)
            c = -65;
            d = 8;
            
            stimThresh=-40%-20;
            fallingThresh=-25%25;
            risingThresh=-20;
            
            vsq=LHANeuron.v(1);
            vl=LHANeuron.v(2);
            cons=LHANeuron.v(3);
     
            p=LHANeuron.fb(1);
            q=LHANeuron.fb(2);
            
            time = LHANeuron.time;
            onset = LHANeuron.onset;
            duration = LHANeuron.duration;
            
            resting = 'resting';
            stimulated = 'stimulated';
            rising = 'rising';
            falling = 'falling';
            
            v=c;
            u= v*p/-q;
            
            Mode=resting
            
            voltages=zeros(time,1);
            currents=zeros(time,1);
            restV=zeros(time,1);
            stimulatedV=zeros(time,1);
            risingV=zeros(time,1);
            fallingV=zeros(time,1);
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
            

            z=('Continuing')
            
            %finding rest A
            
            found=false;
            
            restingIndexs=find(restV);
            stimulatedIndexs=find(stimulatedV);
            risingIndexs=find(risingV)
            fallingIndexs=find(fallingV)
            
            
            vdt=LHANeuron.vdt(2);
            restA=vdt;
            prevRest=0;
            i=1;
            prevStdDif=inf;
            prevAverageDif=inf;
            while(~found)                  
                Ivoltages = zeros(time,1);
                v=c;
                u= v*p/-q;
                
                test=zeros(time,1);

                Mode=resting;
                for t=1:time 
                    I=currents(t);              
                    if(v>=30) 
                       v = c;
                       u = u + d; 
                       
                    end 
                    
                    switch Mode
                       case resting 
                           testRest=~(sum(find(restingIndexs==t))==0);
                           
                           test(t)=testRest;
                                                   
                           if testRest;                               
                              v = v + (v*restA);                             
                              if(I~=0)
                               Mode=stimulated;
                              end                         
                           else
                              v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                              u = u + (p*v + q*u);
                              if(I~=0)
                               Mode=stimulated;
                              end
                           end
                           
                       case stimulated
                            v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                            u = u + (p*v + q*u);
                           if(I==0)&(v<stimThresh)
                               Mode=resting;
                           elseif(v>stimThresh)
                               Mode=rising;
                           end
                       case rising
                            v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                            u = u + (p*v + q*u);
                           if(v>risingThresh)
                               Mode=falling;
                           end
                       case falling
                           v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                           u = u + (p*v + q*u);
                           if(v>=c+fallingThresh)
                               Mode=resting;
                           end
                    end
                    Ivoltages(t)=v;
                end
                i=i+1;
                dif = abs(voltages-Ivoltages);
               
                averageDiff = sum(dif)/time;
                stdDif=std(dif);
                 if(prevAverageDif<=averageDiff)
%                 if(prevStdDif<stdDif)
                    averageDiff=prevAverageDif
                    stdDif=prevStdDif
                    restA=prevRest;
                    found=true;
                 else
                    prevAverageDif=averageDiff;
                    prevStdDif=stdDif;
                    prevRest=restA;
                    restA=restA*.1;
                end
                
            end           
            z='done with rest'
            found=false;
            stimA=vdt;
            prevStim=0;
            i=1;
            prevStdDif=inf;
            prevAverageDif=inf;
            while(~found)                  
                Ivoltages = zeros(time,1);
                v=c;
                u= v*p/-q;
                
                test=zeros(time,1);

                Mode=resting;
                for t=1:time 
                    I=currents(t);              
                    if(v>=30) 
                       v = c;
                       u = u + d;                        
                    end                     
                    switch Mode
                       case resting                           
                          v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                          u = u + (p*v + q*u);            
                          if(I~=0)
                            Mode=stimulated;
                          end                         
                       case stimulated
                           testStim=~(sum(find(stimulatedIndexs==t))==0);
                           test(t)=testStim;
                           if testStim
                               v = v+stimA*v+I;
                               if(I==0)&(v<stimThresh)
                                    Mode=resting;
                               elseif(v>stimThresh)
                                    Mode=rising;
                               end
                           else
                                v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                                u = u + (p*v + q*u);                            
                               if(I==0)&(v<stimThresh)
                                   Mode=resting;
                               elseif(v>stimThresh)
                                   Mode=rising;
                               end
                           end
                       case rising
                            if(v>risingThresh)
                                    v = c;
                                    u = u + d;
                                    Mode=falling;
                            else
                                v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                                u = u + (p*v + q*u);
                           end
                       case falling
                           v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                           u = u + (p*v + q*u);
                           if(v>=c+fallingThresh)
                               Mode=resting;
                           end
                    end
                    Ivoltages(t)=v;
                end
                i=i+1;
                dif = abs(voltages-Ivoltages);
               
                averageDiff = sum(dif)/time;
                stdDif=std(dif);
                 if(prevAverageDif<=averageDiff)
                    averageDiff=prevAverageDif
                    stdDif=prevStdDif;
                    stimA=prevStim;
                    found=true;
                 else
                    prevAverageDif=averageDiff;
                    prevStdDif=stdDif;
                    prevStim=stimA;
                    stimA=stimA*-.1;
                end
                
            end
            z='done with stim'
            i=1;
            prevStdDif=inf;
            prevAverageDif=inf;
            found=false;
            risingA=-vdt;
            prevRising=0;
            while(~found)                  
                Ivoltages = zeros(time,1);
                v=c;
                u= v*p/-q;
                
                test=zeros(time,1);

                Mode=resting;
                for t=1:time 
                    I=currents(t);              
                    if(v>=30) 
                       v = c;
                       u = u + d;                        
                    end                     
                    switch Mode
                       case resting                           
                          v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                          u = u + (p*v + q*u);            
                          if(I~=0)
                            Mode=stimulated;
                          end                         
                       case stimulated        
                                v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                                u = u + (p*v + q*u);                            
                               if(I==0)&(v<stimThresh)
                                   Mode=resting;
                               elseif(v>stimThresh)
                                   Mode=rising;
                               end
                           
                       case rising
                           testRising=~(sum(find(risingIndexs==t))==0);
                           test(t)=testStim;
                           if testRising                              
                              if(v>risingThresh)
                                    v = c;
                                    u = u + d;
                                    Mode=falling;
                               else
                                    v=v+risingA*v+I;  
                               end
                               
                                
                              
                           else
                                if(v>risingThresh)
                                    v = c;
                                    u = u + d;
                                    Mode=falling;
                               else
                                v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                                u = u + (p*v + q*u);
                               end
                           end
                       case falling
                           v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                           u = u + (p*v + q*u);
                           if(v>=c+fallingThresh)
                               Mode=resting;
                           end
                    end
                    Ivoltages(t)=v;
                end
                i=i+1;
                dif = abs(voltages-Ivoltages);
               
                averageDiff = sum(dif)/time;
                stdDif=std(dif);
                 if(prevAverageDif<=averageDiff)
                    averageDiff=prevAverageDif
                    stdDif=prevStdDif
                    risingA=prevRising;
                    found=true;
                 else
                    prevAverageDif=averageDiff
                    prevStdDif=stdDif
                    prevRising=risingA;
                    risingA=abs(prevRising)*.1;
                end
                
            end
            risingA             
            z='done with rising'            
            i=1;
            prevStdDif=inf;
            prevAverageDif=inf;
            found=false;
            fallingA=-vdt;
            prevfalling=0;
            while(~found)                  
                Ivoltages = zeros(time,1);
                v=c;
                u= v*p/-q;
                
                test=zeros(time,1);

                Mode=resting;
                for t=1:time 
                    I=currents(t);              
                    if(v>=30) 
                       v = c;
                       u = u + d;                        
                    end                     
                    switch Mode
                       case resting                           
                          v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                          u = u + (p*v + q*u);            
                          if(I~=0)
                            Mode=stimulated;
                          end                         
                       case stimulated        
                           v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                           u = u + (p*v + q*u);                            
                           if(I==0)&(v<stimThresh)
                               Mode=resting;
                           elseif(v>stimThresh)
                               Mode=rising;
                           end
                           
                       case rising
                            if(v>risingThresh)
                                    v = c;
                                    u = u + d;
                                    Mode=falling;
                            else
                            v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                            u = u + (p*v + q*u);
                            end
                           
                       case falling
                           testFalling=~(sum(find(fallingIndexs==t))==0);
                           test(t)=testFalling;
                           if testFalling
                               v=v+fallingA*v;   
                               if(v>=c+fallingThresh)
                                   Mode=resting;
                               end
                           else
                               v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
                               u = u + (p*v + q*u);
                               if(v>=c+fallingThresh)
                                   Mode=resting;
                               end
                           end
                    end
                    Ivoltages(t)=v;
                end
                i=i+1;
                dif = abs(voltages-Ivoltages);
               
                averageDiff = sum(dif)/time;
                stdDif=std(dif);
                 if(prevAverageDif<=averageDiff)
%                 if(prevStdDif<stdDif)
                    averageDiff=prevAverageDif
                    stdDif=prevStdDif
                    fallingA=prevfalling;
                    found=true;
                 else
                    prevAverageDif=averageDiff;
                    prevStdDif=stdDif;
                    prevfalling=fallingA;
                    fallingA=prevfalling*.1;
                end
                
            end
            z='done with falling'
            z='BEGINING FINAL CHECK'
                
            Ivoltages = zeros(time,1);
            v=c;
            u= v*p/-q;
            test=zeros(time,1);
            Mode=resting
            for t=1:time 
                I=currents(t);              
%                 if(v>=30) 
%                    v = c;
%                    u = u + d;                        
%                 end                     
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
                                        v = c-d;
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

            averageDiff = sum(dif)/time;
            averageError = averageDiff

            figure
            plot(voltages,'b');
            hold on;
            plot(Ivoltages,'r');
            hold on;
            plot(dif,'g');
            hold off;
            xlabel('Time Ms')
            ylabel('Membrane Potential mV')
            axis([0 time -200 200]) 
            
            

           
            
        end
        
        function averageError = generateLHA(LHANeuron)
            c = -65
            d = 8
            
            vsq=LHANeuron.v(1);
            vl=LHANeuron.v(2);
            cons=LHANeuron.v(3);
     
            p=LHANeuron.fb(1);
            q=LHANeuron.fb(2);
            
            time = LHANeuron.time;
            onset = LHANeuron.onset;
            duration = LHANeuron.duration;
            
            firstThres = -1;
            
            v=c;
            u= v*p/-q;
            
            voltages=zeros(time,1);
            I=0; %intial input of zero
            maxV = -inf;
            minV = inf;
            maxVindex = 0;
            minVindex = 0;
            Mode = resting
            for t=1:time %simulate for time
               if t==onset % send dc pulse of I=10 
                   I=10;
               end

               if t==onset+duration
                   I=0;
               end
               if(v>=30) % reset after spike
                   v = c; % reset voltage to rest
                   u = u + d; %reset recovery variable
                   if(firstThres ==-1)
                       firstThres = t 
                   end
               end
               
               v = v + ((vsq*(v^2)) + (vl*v) + cons + -u + I);
               u = u + (p*v + q*u); 
               if v>maxV
                   maxV=v;
                   maxVindex = t;
               end 
               if v<minV
                   minV=v;
                   minVindex=t;
               end
               voltages(t)= v;
            end
            
         
            vdt=LHANeuron.vdt(2)*.001
            
            adj=-70;
            restA=-vdt
            stimA=(voltages(firstThres)-voltages(onset)/(firstThres-onset))*vdt
            risingA=(voltages(maxVindex)-voltages(firstThres)/(maxVindex-firstThres))*vdt*1000
            fallingA=(voltages(minVindex)-voltages(1)/(minVindex-(onset+duration)))*vdt
            
            
            resting = 'resting';
            stimulated = 'stimulated';
            rising = 'rising';
            falling = 'falling';
            
            Ivoltages = zeros(time,1);
            Mode = resting
            I=0; %intial input of zero
            v=c;
            u= v*p/-q;
            for t=1:time %simulate for 1 second
               if t==onset % send dc pulse of I=10 
                   I=10;
               end

               if t==onset+duration
                   I=0;
               end
              
               switch Mode
                   case resting
                       if(I==0)
                            v = ((v + restA*v-u+I));
                            u = u + (p*v + q*u); 
                       else
                           v = ((v + stimA*v -u + I));
                           u = u + (p*v + q*u); 
                           Mode=stimulated
                       end
                   case stimulated
                           v = ((v + stimA*v -u + I));
                           u = u + (p*v + q*u); 
                           if((I==0)&(v<30))
                               Mode=resting
                           elseif(v>=30)
                               Mode=rising
                           end
                   case rising
                            if (v>=maxV)
                                if(I~=0)
                                    v=c;
                                    u = u+d;
                                    Mode=stimulated
                                else
                                    v=minV;
                                    u= u+d;
                                    Mode=falling
                                end
                            else
                                v = ((v + risingA*v -u + I));
                                u = u + (p*v + q*u); 
                            
                            end
                            
                           
                   case falling
                        v = ((v  +fallingA*v)-u+I);
                        u = u + (p*v + q*u); 
                        if(v>=c)
                            Mode=resting
                             t
                             v
                        end
               end
               Ivoltages(t)=v;
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
        end
        
    end
    
end

