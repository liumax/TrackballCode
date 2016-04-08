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
int ttlCount = 0
int intWindow = 200 %this is the window over which counts are integrated

%These are variables for tone triggered TTLs. names describe function.
int pulseNum = 20
int pulseDur = 10
int pulseITI = 40
int pulseCounter = 0

%This is code for cell ID purposes. Toggle holds whether the trip switch
%has been changed, also has TTL pulse information.
int toggle = 0
int idITI = 500
int idPulseDur = 10

clock(slave)

function 1
    while pulseCounter < pulseNum do every pulseITI
        portout[3] = 1
        do in pulseDur
            portout[3] = 0
        end
        pulseCounter = pulseCounter + 1
    then do
        pulseCounter = 0
    end
end;

function 2
    while toggle == 1 do every idITI
        portout[3] = 1
        do in idPulseDur
            portout[3] = 0
    end
end;

callback portin[2] up
	portout[1] = 1
	do in 5
		portout[1] = 0
	end
    ttlCount = ttlCount + 1
    if ttlCount > 1 do
        trigger(1)
    end
    do in intWindow
        ttlCount = ttlCount - 1
    end
end;

callback portin[4] up
    if toggle == 0 do
        toggle = 1
        trigger(2)
    elseif toggle == 1 do
        toggle = 0
    end
end;



