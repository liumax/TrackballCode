%This is code to fix problems where TTLs are not delivered/detected. This
%is meant to use the predefined ITI times and number of expected inputs to
%fix this problem. This generates a dummy time that is used to keep code
%functional. Also plants a warning that indicates that there was a missed
%pulse. 


function [s,repairedTTLs] = functionTTLRepairSystem(numExpected,expectedITIs,actualTimes,pairingToggle,laserSig,laserLag,s);

%% Parameters:
defaultLag = 0.25; %default lag in seconds


%% Compute differences from actual TTL inputs
actualDiffs = diff(actualTimes);

%% Split analysis based on pairing vs non-pairing:
if pairingToggle == 1
%     try
%         laserSig = soundFile.Pairing.LaserTriggerPulseITI;
%     catch
%         disp('Cant Find Laser Pulse ITI, Defaulting to 4ms')
%         laserSig = 0.004; %hardcoded fix in case this doesnt pull things up
%     end
%     try
%         laserLag = soundFile.Pairing.OptoStimDelay;
%     catch
%         disp('Cant Find Laser Delay, Defaulting to -4ms')
%         laserLag = -0.004; %hardcoded fix if cant pull up appropriate value
%     end
    
    %find short spacing. This should be either lasers, or target.
    shortFinder = find(actualDiffs <= laserSig*2);

    %now I need to try and remove the excess TTLs that are for the laser.
    if laserSig == -laserLag;
        %in this case, there should be two ttls for each paired tone and one
        %for each control. Therefore, the total number should be 3x the
        %expected number
        totalExpected = 3*numExpected;
        %next, identify and separate out the TTLs based on their spacing.
        %Target TTLs should have the short spacing, controls have longer
        %spacing
        laserTTLs = actualTimes(shortFinder);
        %eliminate laser TTLs from actualTimes
        remainingTTLs = actualTimes;
        remainingTTLs(shortFinder) = [];
    elseif laserSig ~= -laserLag;
        %in this case, laser is completely separate from the tone. As a result,
        %there should be 4x the expected number
        totalExpected = 4*numExpected;

        %identify and separate out TTLs based on spacing
        laserTTLs = actualTimes(sort([shortFinder;shortFinder+1]));
        %eliminate laser TTLs (remove both)
        remainingTTLs = actualTimes;
        remainingTTLs(sort([shortFinder;shortFinder+1])) = [];
    end
else
    remainingTTLs = actualTimes;
end

%% Recompute differences
remainingDiff = diff(remainingTTLs);

%% Run while loop to execute repair code. 
checkLoop = 0;
startInd = 1;
repairInd = 1;
while checkLoop == 0;
    %check to see if values are roughly correct
    try
        valCheck = remainingDiff(startInd) - expectedITIs(startInd);
        if abs(valCheck) <= defaultLag %this means that the value is within expectations
            startInd = startInd + 1;
        else %this means that the ITI isnt good looking! XD
            disp('Aberrant ITI Time Detected')
            %determine if ITI time is too long or short
            if valCheck > defaultLag %too long
                disp('TTL Missing')
                %now I want to add a value based on estimated value
                remainingTTLs(startInd + 2:end+1) = remainingTTLs(startInd+1:end);
                %compute average delay
                delayCalc = mean(remainingDiff(1:startInd)-expectedITIs(1:startInd));
                remainingTTLs(startInd+1) = remainingTTLs(startInd) + expectedITIs(startInd) + delayCalc;
                remainingDiff = diff(remainingTTLs);
                disp('TTL Replaced Based on ITI')
                s.TTLRepair.(targetName).Flag{repairInd} = 'MissingTTL';
                s.TTLRepair.(targetName).Index{repairInd} = startInd;
                repairInd = repairInd + 1;
                %dont update start ind. want to see if things were fixed.
            elseif valCheck < defaultLag; %too short
                %remove the extraneous TTL
                disp('Extraneous TTL Detected, Removing')
                remainingTTLs(startInd+1) = [];
                remainingDiff = diff(remainingTTLs);
                s.TTLRepair.(targetName).Flag{repairInd} = 'ExtraTTL';
                s.TTLRepair.(targetName).Index{repairInd} = startInd;
                repairInd = repairInd + 1;
                %dont update start ind, since I just deleted a value
            end
        end
    catch
        %first check if things are working!
        lengthCheck = length(remainingTTLs) - length(expectedITIs);
        if lengthCheck == 0
            disp('Damage Repaired')
            checkLoop = 1;
        else
            error('Error in While Loop')
        end
    end
end

repairedTTLs = remainingTTLs;



end


















