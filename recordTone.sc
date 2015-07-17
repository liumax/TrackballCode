%This is the sc file
%ports 1 is output to sound
%port 2 is output to plexon (habituation)
%port 3 is output to plexon (recording)
%port 4 is laser trigger
%port 5 is actual laser output

%port 8 is the trigger to go from normalizing trials (getting stable baseline)
%to recording trials.

int toneDur
int stimDur
int stimFreq
int stimNum
int stimTime
int stimCounter = 0

int itiDur

callback portin[4] up
    portout[4] = 0
    if stimCounter < stimNum do
        portout[5] = 1
        do in stimDur
            portout[5] = 0
            stimCounter = stimCounter + 1
        end
        do in stimFreq
            portout[4] = 1
        end
    else do
        stimCounter = 0
        portout[4] = 0
    end
end;

callback portin[8] up
    disp('Habituation Trials Ended')
    disp('Recording Trials Begin Now')
end;

function 1
    do in itiDur
        disp('Sound Start')
        portout[1] = 1
        portout[2] = 1
        do in toneDur
            portout[1] = 0
            portout[2] = 0
            disp('Sound End')
            disp('TriggerMatlab')
        end
    end
end;

function 2
    do in itiDur
        disp('Laser Stim!')
        portout[4] = 1
        do in stimTime
            disp('Sound Start')
            portout[1] = 1
            portout[3] = 1
            do in toneDur
                portout[1] = 0
                portout[3] = 0
                disp('Sound End')
                disp('TriggerMatlab')
            end
        end
    end
end;



