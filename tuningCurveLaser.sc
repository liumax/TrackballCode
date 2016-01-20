%this is code for reading out the TTLs from the Roland
%The goal is to detect all TTLs and respond to a subset of them.
%Currently, it reads from the TTLs, and responds if a double TTL signal occurs

int ttlCount = 0
int intWindow = 1000 %this is the window over which counts are integrated
int laserNum 
int laserDur 
int laserITI 
int pulseCounter = 0

%function 1 generates laser pulse
function 1
    while pulseCounter < laserNum do every laserITI
        portout[2] = 1
        do in laserDur
            portout[2] = 0
        end
        pulseCounter = pulseCounter + 1
    then do
        pulseCounter = 0
    end
end;

callback portin[1] up
    ttlCount = ttlCount + 1
    if ttlCount > 1 do
        trigger(1)
    end
    do in intWindow
        ttlCount = ttlCount - 1
    end
end;





