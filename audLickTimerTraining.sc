%This is the sc file
%port 1 is lickometer
%port 2 is reward solenoid
%ports 3-4 are sounds. 3 is low tone, 4 is high tone. 
%port 5 is warning light telling mouse they are licking in between

int cueID %determines whether to play high or low tones
int cueDur %determines how long in ms cue lasts
int itiDur %duration of intertrial interval
int rewDur %reward duration
int rewWin %reward window (should be 2 seconds)
int lickTicker = 0 %this is self-updating counter of licks
int noLick %this is period over which no licking is allowed

int rewSwitch = 0

callback portin[1] up
    disp('Lick Detected')
    lickTicker = lickTicker + 1
    do in noLick
        lickTicker = lickTicker - 1
    end
    if rewSwitch == 1 do
        rewSwitch = 0
        disp('Reward Delivered')
        portout[2] = 1
        do in rewDur
            portout[2] = 0
            do in itiDur
                disp('TriggerMatlab')
            end
        end
    end
end;

function 2
    do in noLick
        if lickTicker > 0 do
            disp('Bad Licks')
            portout[5] = 1
            do in rewDur
                portout[5] = 0
            end
            trigger(2)
        else if lickTicker == 0 do
            disp('Initiating Reward Trial')
            disp('Cue On')
            portout[cueID] = 1
            rewSwitch = 1
            do in cueDur
                portout[cueID] = 0
            end
            do in rewWin
                if rewSwitch == 1 do
                    rewSwitch = 0
                    do in itiDur
                        disp('TriggerMatlab')
                    end
                end
            end
        end
    end
end;


function 1
    if lickTicker > 0 do
        disp('Bad Licks')
        portout[5] = 1
        do in rewDur
            portout[5] = 0
        end
        trigger(2)
    else if lickTicker == 0 do
        disp('Initiating Reward Trial')
        disp('Cue On')
        portout[cueID] = 1
        rewSwitch = 1
        do in cueDur
            portout[cueID] = 0
        end
        do in rewWin
            if rewSwitch == 1 do
                rewSwitch = 0
                do in itiDur
                    disp(itiDur)
                    disp('TriggerMatlab')
                end
            end
        end
    end
end;



