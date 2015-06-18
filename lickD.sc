%This is the sc file
%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int rewLength % CHANGE THIS IN BLOCKS
int soundRewDel % CHANGE THIS EVERY TRIAL
int soundDur %maybe hardcode this? doesnt change.
int preDelay = 1000
int postDelay = 2000
int trialState
int itiDur

%These are variables for tracking running disk.
int upA = 0
int velCount = 0
int intWindow1 = 500

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
                disp(rewLength)
                portout[4] = 1
                    do in rewLength
                        portout[4] = 0
                        disp('Reward Completed')
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

callback portin[7] up 
    upA = upA + 1
    if upA == 3 do
        velCount = velCount + 1
        disp(velCount)
        upA = 0
        do in intWindow1
            velCount = velCount - 1
            disp(velCount)
        end
    end
end;