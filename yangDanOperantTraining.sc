%This is the sc file
%port 1 is lickometer
%port 2 is reward solenoid
%port 3 is punishment solenoid 
%ports 4-6 are lights. 4 warning 5 reward 6 punishment

int warning
int warningDelay
int cueDur
int graceDur
int itiDur
int rewDur

int rewSwitch = 0


callback portin[1] up
    disp('Lick Detected')
    if rewSwitch == 1 do
        rewSwitch = 0
        disp('Reward Delivered')
        portout[2] = 1
        do in rewDur
            portout[2] = 0
            do in warningDelay
                portout[5] = 0
            end
            do in itiDur
                disp('TriggerMatlab')
            end
        end
    end
end


function 1
    disp('Initiating Reward Trial')
    disp('Warning Light On')
    portout[4] = 1
    do in warning
        disp('Warning Light Off')
        portout[4] = 0
    end
    do in warningDelay
        disp('Cue Light On')
        portout[5] = 1
        rewSwitch = 1
    end
end;



