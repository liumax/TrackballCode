%This is the sc file for pavTwoToneCatchFreeScript, which is meant to be a 
%script to perform two tone pavlov behavior along with free rewards AND catch trials. 

%This has only one sound (port 2)
%This has lickometer at port 3, solenoid at port 4

int lickWind %This is for licking isolation window. Not used here, but need variable so things dont break
int rewLength %duration of reward, will change every trial
int toneRewDel % delay between tone and reward. fixed for entire session
int signalDel %delay between reward and triggering matlab again. 

int timeDelay %delay to trigger matlab callback later.
int itiDur
int trialMode = 0 %is something to lock things so that extraneous signals dont fuck shit up

int laserDur = 0 %initializes laser duration at zero so doesnt throw bug


function 1 %This function merely serves to wait out the ITI
    disp('TrialStart')
    do in itiDur
	trialMode = 1
	disp('TriggerSound')
    end
end;

function 2 %This function will be for triggering the laser
    disp('LaserTrial')
    portout[6] = 1
    do in laserDur
        portout[6] = 0
	disp('Laser Delivered')
	do in toneRewDel
		do in rewLength
            do in 5000
                disp('PlotTime')
                do in 1000
                    disp('TriggerMatlab')
                end
            end
		end
	end
    end
end;

callback portin[1] up %tone has been detected
	portout[8] = 1
	do in 5
		portout[8] = 0
	end
    disp('Tone Delivered')
	if trialMode == 1 do in 0
		trialMode = 0
		do in toneRewDel
			portout[2] = 1 %deliver reward
			do in rewLength
				portout[2] = 0
                do in 5000
                    disp('PlotTime')
                    do in 1000
                        disp('TriggerMatlab') %start sequence over again!
                    end
                end
			end
		end
	end
end;

callback portin[2] up %lick has been detected
    disp('Lick Detected')
end;

