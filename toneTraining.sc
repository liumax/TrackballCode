%This is the sc file
%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int rewLength 
int soundRewDel % CHANGE THIS EVERY TRIAL
int soundDur 
int sound

int lickWindow %another delay, just convenient.
int timeDelay %delay to trigger matlab callback later.
int itiDur

function 1
    disp('Initiating trial')
    do in itiDur
        do in lickWindow
            if sound == 0 do
                portout[1] = 1 % sound on, 
                disp('SoundOn')
                disp('White Noise')
                do in soundDur
                    portout[1] = 0 % sound off
                    disp('SoundOff')
                    do in timeDelay
                        disp('TriggerMatlab')
                    end
                end
            else if sound == 1 do
                portout[2] = 1 % sound on
                disp('SoundOn')
                disp('Pure Tone')
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
            end
        end
    end
end;

function 3
    disp('Initiating trial')
    do in itiDur
        do in lickWindow
            if sound == 0 do
                portout[2] = 1 % sound on, 
                disp('SoundOn')
                disp('Pure Tone')
                do in soundDur
                    portout[2] = 0 % sound off
                    disp('SoundOff')
                    do in timeDelay
                        disp('TriggerMatlab')
                    end
                end
            else if sound == 1 do
                portout[1] = 1 % sound on
                disp('SoundOn')
                disp('White Noise')
                do in soundDur
                    portout[1] = 0 % sound off
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
            end
        end
    end
end;

