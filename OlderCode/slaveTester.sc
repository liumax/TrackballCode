%this is code for reading out the TTLs from the Roland
%The goal is to detect all TTLs and respond to a subset of them.

int ttlCount = 0
int intWindow = 1000 %this is the window over which counts are integrated
int pulseNum = 20
int pulseDur = 10
int pulseITI = 40
int pulseCounter = 0

clock(slave)

function 1
    while pulseCounter < pulseNum do every pulseITI
        portout[2] = 1
        do in pulseDur
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





