%This is the sc file for twoTonePavlovWAV. This is meant to run the behavior with a WAV trigger
%Port 2: solenoid and lick detection
%Port 1: sound! goes to wav trigger

%NEED TO REASSIGN PORTS BASED ON WHAT IS AVAILABLE ON MBED
%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int lickWind %This is for licking isolation window. Not used here, but need variable so things dont break
int rewLength %duration of reward, will change every trial
int toneRewDel % delay between tone and reward. fixed for entire session
int signalDel %delay between reward and triggering matlab again. 

int timeDelay %delay to trigger matlab callback later.
int itiDur

int laserDur = 0 %initializes laser duration at zero so doesnt throw bug

int trPhase = 0%0 means pre tone delivery
                %1 means during tone
                %2 means 5 seconds post tone
                %4 means during laser
                %5 means post laser (5sec)

function 1 %This function waits out ITI and then delivers sound and reward
    disp('SoundTrial')
    do in itiDur
        trPhase = 1
        disp(trPhase)
        portout[1] = 1
        do in
            portout[1] = 0
        end
        do in toneRewDel
            portout[2] = 1
            do in rewLength
                portout[2] = 0
                trPhase = 2
                disp(trPhase)
                do in 5000
                    trPhase = 0
                    disp(trPhase)
                    disp('PlotTime')
                    do in 1000
                        disp('TriggerMatlab')
                    end
                end
            end
        end
    end
end;

function 2 %This function will be for triggering the laser
    disp('LaserTrial')
    trPhase = 1
    disp(trPhase)
    portout[6] = 1
    do in laserDur
        trPhase = 4
        disp(trPhase)
        portout[6] = 0
	disp('Laser Delivered')
	do in toneRewDel
		do in rewLength
            trPhase = 5
            disp(trPhase)
            do in 5000
                trPhase = 0
                disp(trPhase)
                disp('PlotTime')
                do in 1000
                    disp('TriggerMatlab')
                end
            end
		end
	end
    end
end;

callback portin[2] up %lick has been detected
    disp('Lick Detected')
end;

