%This is the sc file
%This has two sounds (ports 1&2)
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
int sound

function 1
    disp('Initiating trial')
    trialState = 4
    do in itiDur
        trialState = 1
        do in preDelay
            trialState = 2
            if sound == 1
                portout[1] = 1 % sound on
            elseif sound == 2
                portout[2] = 1 % sound on
            end
            disp('SoundOn')
            do in soundDur
                if sound == 1
                    portout[1] = 0 % sound off
                elseif sound == 2
                    portout[2] = 0 % sound off
                end
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