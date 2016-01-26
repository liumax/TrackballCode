%This is the sc file
%ports 1 and 2 are inputs from other MBED
%port 3 is solenoid
%port 4 is LED to signal stop

int rewLength

int timeDelay %delay to trigger matlab callback later.
int itiDur
int waitPeriod

int motion = 0
int stopped = 0

int moveSwitch = 0
int trialSwitch = 0

callback portin[1] up
    motion = 1
    if trialSwitch == 1 do
        trialSwitch = 0
        moveSwitch = 1
        disp('TriggerMatlab')
    end
end;

callback portin[2] up
    stopped = 1
    disp('Stopped')
end;

callback portin[1] down
    motion = 0
end;

callback portin[2] down
    stopped = 0
end;

function 2
    disp('Failed to Stop')
    do in waitPeriod
        if motion == 0 && stopped == 1 do
            disp('Successful Stop')
            disp('Reward Delivered')
            portout[3] = 1
            do in rewLength
                portout[3] = 0
            end
            portout[4] = 0
            do in itiDur
                moveSwitch = 0
                trialSwitch = 1
            end
        else do
            trigger(2)
        end
    end
end;


function 1
    disp('Initiating trial')
    disp('Light On')
    portout[4] = 1
    do in waitPeriod
        if motion == 0 && stopped == 1 do
            disp('Successful Stop')
            disp('Reward Delivered')
            portout[3] = 1
            do in rewLength
                portout[3] = 0
            end
            portout[4] = 0
            do in itiDur
                moveSwitch = 0
                trialSwitch = 1
            end
        else do
            trigger(2)
        end
    end
end;


function 3
    disp('Initiating trial')
    do in itiDur
        disp('Light On')
        portout[4] = 1
        do in waitPeriod
            if motion == 0 && stopped == 1 do
                disp('Successful Stop')
                disp('Reward Delivered')
                portout[3] = 1
                do in rewLength
                    portout[3] = 0
                end
                portout[4] = 0
                do in itiDur
                    moveSwitch = 0
                    trialSwitch = 1
                end
            else do
                trigger(2)
            end
        end
    end
end;
