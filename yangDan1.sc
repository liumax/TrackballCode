%This is the sc file
%port 1 is lickometer
%port 2 is reward solenoid
%port 3 is punishment solenoid 
%ports 4-7 are lights

int warning
int warningDelay
int cueDur
int graceDur
int itiDur
int rewDur
int punDur

int rewSwitch = 0
int punSwitch = 0

callback portin[1] up
    disp('Lick Detected')
    if rewSwitch == 1 do
        disp('Reward Delivered')
        rewSwitch = 0
        portout[2] = 1
        do in rewDur
            disp('Reward Terminated')
            portout[2] = 0
        end
    end
    if punSwitch == 1 do
        disp('Punishment Delivered')
        punSwitch = 0
        portout[3] = 1
        do in punDur
            disp('Punishment Terminated')
            portout[3] = 0
        end
    end
end;

function 1
    disp('Initiating Neutral Trial')
    disp('Warning Light On')
    portout[4] = 1
    do in warning
        disp('Warning Light Off')
        portout[4] = 0
    end
    do in warningDelay
        disp('Cue Light On')
        portout[5] = 1
        do in graceDur
            disp('Grace Period Ended')
        end
        do in cueDur
            disp('Cue Light Off')
            portout[5] = 0
            do in itiDur
                disp('TriggerMatlab')
            end
        end
    end
end;

function 2
    disp('Initiating Reward Trial')
    disp('Warning Light On')
    portout[4] = 1
    do in warning
        disp('Warning Light Off')
        portout[4] = 0
    end
    do in warningDelay
        disp('Cue Light On')
        portout[6] = 1
        do in graceDur
            disp('Grace Period Ended')
            rewSwitch = 1
        end
        do in cueDur
            disp('Cue Light Off')
            portout[6] = 0
            rewSwitch = 0
            do in itiDur
                disp('TriggerMatlab')
            end
        end

    end
end;

function 3
    disp('Initiating Punishment Trial')
    disp('Warning Light On')
    portout[4] = 1
    do in warning
        disp('Warning Light Off')
        portout[4] = 0
    end
    do in warningDelay
        disp('Cue Light On')
        portout[7] = 1
        do in graceDur
            disp('Grace Period Ended')
            punSwitch = 1
        end
        do in cueDur
            disp('Cue Light Off')
            portout[7] = 0
            punSwitch = 0
            do in itiDur
                disp('TriggerMatlab')
            end
        end
    end
end;


