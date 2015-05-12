%This is the sc file
%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int rewLength % CHANGE THIS IN BLOCKS
int soundRewDel % CHANGE THIS EVERY TRIAL
int baitDur
int bait1 = 0
int soundDur
int preDelay = 1000
int postDelay = 2000
int trialState
int itiDur

function 1
    disp('Initiating trial')
    trialState = 4
    do in itiDur
        trialState = 1
        do in preDelay
            trialState = 2
            portout[2] = 1 % sound on
            disp('SoundOn')
            do in soundDur
                portout[2] = 0 % sound off
                disp('SoundOff')
            end
            do in soundRewDel
                trialState = 3
                disp('Reward Delivered')
                portout[4] = 1
                    do in rewLength
                        portout[4] = 0
                        disp(rewLength)
                        do in postDelay
                            trialState = 4
                        end
                    end
            end
        end
    end
end;

callback portin[3] up % lickometer activated
    disp(trialState)
end;