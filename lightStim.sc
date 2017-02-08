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
int idITI = 2000
int idPulseDur = 500
int idCounter = 0

function 1
	while toggle == 1 do every idITI
		portout[6] = 1
		do in idPulseDur
			portout[6] = 0
		end
		idCounter = idCounter + 1
		disp(idCounter)
	end
end;


callback portin[4] up
    if toggle == 0 do
        toggle = 1
        trigger(1)
    else do
        toggle = 0
	idCounter = 0
	disp('Stop ID')
    end
end;



