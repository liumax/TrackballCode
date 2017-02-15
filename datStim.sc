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

int ITI = 0;

%These are variables for tone triggered TTLs. names describe function.
%These variables are for dopamine pulsing
int pulseNum = 20
int pulseDur = 10
int pulseITI = 30
int pulseCounter = 0


function 1
    do in ITI
        portout[1] = 1 %signals start of pulse
        do in pulseDur
            portout[1] = 0
        end
        while pulseCounter < pulseNum do every pulseITI
            portout[5] = 1
            do in pulseDur
                portout[5] = 0
            end
            pulseCounter = pulseCounter + 1
        then do
            pulseCounter = 0
            disp('TriggerMatlab')
        end
    end
end;



