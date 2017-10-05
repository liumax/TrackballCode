%This code is meant to solely execute manual stimulation, with no passing of auditory TTLs or other signals. 

%These are variables for tone triggered TTLs. names describe function.
%These variables are for dopamine pulsing
int pulseNum = 100
int pulseDur = 2000
int pulseITI = 6000
int pulseCounter = 0

%This is code for cell ID purposes. Toggle holds whether the trip switch
%has been changed, also has TTL pulse information.
int toggle = 0


function 1
	portout[1] = 1
	do in pulseDur
		portout[1] = 0
	end
   	while pulseCounter < pulseNum do every pulseITI
        	portout[6] = 1
       	 do in pulseDur
            		portout[6] = 0
        	end
       	 pulseCounter = pulseCounter + 1
		disp(pulseCounter)
   	 then do
        	pulseCounter = 0
		disp('Stimulation Finished')
    	 end
end;

callback portin[3] up
    if toggle == 0 do
        toggle = 1
        trigger(1)
    else do
        toggle = 0
	disp('Stop Stim')
    end
end;
