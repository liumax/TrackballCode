%This is the sc file
%port 1 is to MCU
%port 2 is input from Roland
%port 3 is output to laser
%port 4 is output to wavTrigger. 

%port 8 is the trigger to go from normalizing trials (getting stable baseline)
%to recording trials.

int trigDur

int itiDur

%function 1 generates noise!
function 1
    do in itiDur
        disp('Sound Start')
        portout[1] = 1
        portout[4] = 1
        do in trigDur
            portout[1] = 0
            portout[4] = 0
            disp('Sound End')
            disp('TriggerMatlab')
        end
    end
end;




