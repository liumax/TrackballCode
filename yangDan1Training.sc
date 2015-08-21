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

int rewSwitch = 0
int punSwitch = 0

callback portin[1] up
    disp('Lick Detected')
end

function 1
    disp('Initiating Training Trial')
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
            disp('Reward Delivered')
            portout[2] = 1
            do in rewDur
                portout[2] = 0
            end
        end
        do in cueDur
            disp('Cue Light Off')
            portout[6] = 0
            do in itiDur
                disp('TriggerMatlab')
            end
        end
    end
end;


