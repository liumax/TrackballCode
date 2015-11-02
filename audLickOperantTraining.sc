%This is the sc file
%port 1 is lickometer
%port 2 is reward solenoid
%ports 3-4 are sounds. 3 is low tone, 4 is high tone. 

int cueID %determines whether to play high or low tones
int cueDelay %determines how long in ms cue lasts past reward delivery
int itiDur %duration of intertrial interval
int rewDur %reward duration

int rewSwitch = 0

callback portin[1] up
    disp('Lick Detected')
    if rewSwitch == 1 do
        rewSwitch = 0
        disp('Reward Delivered')
        portout[2] = 1
        do in rewDur
            portout[2] = 0
            do in cueDelay
                portout[cueID] = 0
            end
            do in itiDur
                disp('TriggerMatlab')
            end
        end
    end
end


function 1
    disp('Initiating Reward Trial')
    disp('Cue Light On')
    portout[cueID] = 1
    rewSwitch = 1
end;

function 2
    disp('Initiating Reward Trial')
    do in itiDur
        disp('Cue Light On')
        portout[cueID] = 1
        rewSwitch = 1
    end
end;



