%this is code for reading out the TTLs from the Roland
%The goal is to detect all TTLs and respond to a subset of them.
%output 1 is devoted to sending TTL pulses to MCU regarding tone timing
%input 2 is reading out the Roland outputs
%output 3 is devoted to laser output!
%input 4 is going to be a manual toggle switch! This will change between
%modes of passive function (listening and passing on tone TTLs) and active
%(Cell ID protocol).

%This is code for integrating double pulses for tone triggered lasers. TTL
%count is a holder for leaky integration. intWindow is the time over which
%leaky integrator works. 

%This is code for cell ID purposes. Toggle holds whether the trip switch
%has been changed, also has TTL pulse information.
int toggle = 0
int ITI = 0;
int idPulseDur = 500
int idCounter = 0

%These are variables for tone triggered TTLs. names describe function.
%These variables are for dopamine pulsing
int pulseNum = 20
int pulseDur = 10
int pulseITI = 30
int pulseCounter = 0


function 1
    while pulseCounter < pulseNum do every pulseITI
        portout[5] = 1
        do in pulseDur
            portout[5] = 0
        end
        pulseCounter = pulseCounter + 1
    then do
        pulseCounter = 0
    end
end;

function 2
    if toggle == 1
        ITI = ITImin+
        do in 
        end
    end
end;

callback portin[4] up
    if toggle == 0 do
        toggle = 1
        trigger(2)
    else do
        toggle = 0
	idCounter = 0
	disp('Stop ID')
    end
end;



