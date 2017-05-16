%This is the sc file for twoTonePavlovComputerTraining. This is meant to be
%performed with the spout in the mouse's mouth. This means that you basically
%are not collecting any licking data. 

%NEED TO REASSIGN PORTS BASED ON WHAT IS AVAILABLE ON MBED
%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int lickWind %This is for licking isolation window. Not used here, but need variable so things dont break
int rewLength %duration of reward, will change every trial
int toneRewDel % delay between tone and reward. fixed for entire session
int signalDel %delay between reward and triggering matlab again. 

int timeDelay %delay to trigger matlab callback later.
int itiDur

function 1 %This function merely serves to wait out the ITI
    disp('Initiating trial')
    do in itiDur
        disp('TriggerSound')
    end
end;

callback portin[] up %tone has been detected
    disp('Tone Delivered')
    do in toneRewDel
        portout[] = 1 %deliver reward
        do in rewLength
            portout[] = 0
        end
        do in 
            disp('TriggerMatlab') %start sequence over again!
        end
    end
end;
