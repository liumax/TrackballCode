%This is the sc file
%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int rewLength % CHANGE THIS IN BLOCKS
int soundRewDel % CHANGE THIS EVERY TRIAL
int soundDur %maybe hardcode this? doesnt change.

int timeDelay %delay to trigger matlab callback later.
int itiDur

int lickCounter = 0
int lickWindow


function 2
    disp('Bad Licks!')
    do in lickWindow
        if lickCounter < 1 do
            portout[2] = 1 % sound on
            disp('SoundOn')
            do in soundDur
                portout[2] = 0 % sound off
                disp('SoundOff')
                do in timeDelay
                    disp('TriggerMatlab')
                end
            end
            do in soundRewDel
                if lickCounter > 0 do
                    disp('Reward Delivered')
                    disp(rewLength)
                    portout[4] = 1
                        do in rewLength
                            portout[4] = 0
                            disp('Reward Completed')
                        end
                else do
                    disp('No Reward')
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
        do in lickWindow
            if lickCounter < 1 do
                portout[2] = 1 % sound on
                disp('SoundOn')
                do in soundDur
                    portout[2] = 0 % sound off
                    disp('SoundOff')
                    do in timeDelay
                        disp('TriggerMatlab')
                    end
                end
                do in soundRewDel
                    if lickCounter > 0 do
                        disp('Reward Delivered')
                        disp(rewLength)
                        portout[4] = 1
                            do in rewLength
                                portout[4] = 0
                                disp('Reward Completed')
                            end
                    else do
                        disp('No Reward')
                    end
                end
            else do
                trigger(2)
            end
        end
    end
end;

function 3 %This is for training trials.
    portout[2] = 1 % sound on
    disp('SoundOn')
    do in soundDur
        portout[2] = 0 % sound off
        disp('SoundOff')
        do in timeDelay
            disp('TriggerMatlab')
        end
    end
    do in soundRewDel
        disp('Reward Delivered')
        disp(rewLength)
        portout[4] = 1
            do in rewLength
                portout[4] = 0
                disp('Reward Completed')
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

