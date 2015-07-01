%This is the sc file
%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int rewLength % CHANGE THIS IN BLOCKS
int soundRewDel % CHANGE THIS EVERY TRIAL
int soundDur %maybe hardcode this? doesnt change.

int lickWindow

int timeDelay %delay to trigger matlab callback later.
int itiDur

function 1
    disp('Initiating trial')
    do in itiDur
        do in lickWindow
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
        end
    end
end;