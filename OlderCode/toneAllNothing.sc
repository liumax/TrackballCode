%This is the sc file
%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int rewLength 
int soundRewDel % CHANGE THIS EVERY TRIAL
int soundDur 
int sound
int toneRule

int timeDelay %delay to trigger matlab callback later.
int itiDur

int lickCounter = 0
int lickWindow


function 2
    disp('Bad Licks!')
    do in lickWindow
        if lickCounter < 1 do
            if sound == 0 do
                portout[1] = 1 % sound on
                disp('SoundOn')
                disp('White Noise')
            else if sound == 1 do
                portout[2] = 1 % sound on
                disp('SoundOn')
                disp('Pure Tone')
            end
            do in soundDur
                if sound == 0 do
                    portout[1] = 0 % sound off
                else if sound == 1 do
                    portout[2] = 1 % sound on
                end
                disp('SoundOff')
                do in timeDelay
                    disp('TriggerMatlab')
                end
            end
            do in soundRewDel
                if sound == 0 && toneRule == 0 do
                    if lickCounter > 0 do
                        disp('Reward Delivered')
                        portout[4] = 1
                        do in rewLength
                            portout[4] = 0
                            disp('Reward Completed')
                        end
                    else do
                        disp('No Reward')
                    end
                else if sound == 1 && toneRule == 1 do
                    if lickCounter > 0 do
                        disp('Reward Delivered')
                        portout[4] = 1
                        do in rewLength
                            portout[4] = 0
                            disp('Reward Completed')
                        end
                    else do
                        disp('No Reward')
                    end
                else do
                    disp('Unrewarded')
                end
            end
        else do
            trigger(2)
        end
    end
end;




function 1
    disp('Initiating trial')
    do in itiDur
        if lickCounter < 1 do
            if sound == 0 do
                portout[1] = 1 % sound on
                disp('SoundOn')
                disp('White Noise')
            else if sound == 1 do
                portout[2] = 1 % sound on
                disp('SoundOn')
                disp('Pure Tone')
            end
            do in soundDur
                if sound == 0 do
                    portout[1] = 0 % sound off
                else if sound == 1 do
                    portout[2] = 1 % sound on
                end
                disp('SoundOff')
                do in timeDelay
                    disp('TriggerMatlab')
                end
            end
            do in soundRewDel
                if sound == 0 && toneRule == 0 do
                    if lickCounter > 0 do
                        disp('Reward Delivered')
                        portout[4] = 1
                        do in rewLength
                            portout[4] = 0
                            disp('Reward Completed')
                        end
                    else do
                        disp('No Reward')
                    end
                else if sound == 1 && toneRule == 1 do
                    if lickCounter > 0 do
                        disp('Reward Delivered')
                        portout[4] = 1
                        do in rewLength
                            portout[4] = 0
                            disp('Reward Completed')
                        end
                    else do
                        disp('No Reward')
                    end
                else do
                    disp('Unrewarded')
                end
            end
        else do
            trigger(2)
        end
    end
end;


callback portin[3] up % lickometer activated
    lickCounter = lickCounter + 1
    disp('Lick Detected')
    do in lickWindow
        lickCounter = lickCounter - 1
    end
end;

