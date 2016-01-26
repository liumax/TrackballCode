%This is the sc file
%port 1 is lickometer
%port 2 is reward solenoid
%ports 3-4 are sounds. 3 is low tone, 4 is high tone. 
%port 5 is warning light telling mouse they are licking in between

int cueID %determines whether to play high or low tones
int distID %determines the distractor tone (opposite of cueID)
int cueDur %determines how long in ms cue lasts
int itiDur %duration of intertrial interval
int rewDur %reward duration
int rewWin %reward window (should be 2 seconds)
int lickTicker = 0 %this is self-updating counter of licks
int noLick %this is period over which no licking is allowed
int trialType %This instructs which trial type to execute

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

function 3
    disp('Initiating Reward Trial')
    disp('Reward Cue On')
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
end;

function 4
    disp('Initiating Distractor Trial')
    disp('Distractor Cue On')
    portout[distID] = 1
    do in cueDur
        portout[distID] = 0
        do in itiDur
            disp(itiDur)
            disp('TriggerMatlab')
        end
    end
end;

function 5
    disp('Initiating Uncued Reward')
    disp('Uncued Reward Delivered')
    portout[2] = 1
    do in rewDur
        portout[2] = 0
        do in itiDur
            disp(itiDur)
            disp('TriggerMatlab')
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
            if trialType == 3 do
                trigger(3)
            else if trialType == 4 do
                trigger(4)
            else if trialType == 5 do
                trigger(5)
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
        if trialType == 3 do
            trigger(3)
        else if trialType == 4 do
            trigger(4)
        else if trialType == 5 do
            trigger(5)
        end
    end
end;

function 6 %this is here for the first trial so there is a delay
    do in itiDur
        if lickTicker > 0 do
            disp('Bad Licks')
            portout[5] = 1
            do in rewDur
                portout[5] = 0
            end
            trigger(2)
        else if lickTicker == 0 do
            if trialType == 3 do
                trigger(3)
            else if trialType == 4 do
                trigger(4)
            else if trialType == 5 do
                trigger(5)
            end
        end
    end
end;




