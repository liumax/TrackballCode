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
    end
end;

callback portin[1] up %tone has been detected
	if trialMode == 1 do in 0
		trialMode = 0
		disp('Tone Delivered')
		do in toneRewDel
			portout[2] = 1 %deliver reward
			do in rewLength
				portout[2] = 0
		%disp(toneRewDel)
		%disp(rewLength)
				disp('TriggerMatlab') %start sequence over again!
			end
		end
	end
end;
