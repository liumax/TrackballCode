%this is meant to be a function that will help to align TDT and MBED pulse
%times and fix them if there are errors. Specifically, this is designed for
%the jitter signal, which will often have many more events in the MBED file
%as opposed to the TDT file, since it will also include samples from the
%testing period

function [tdtTimes,mbedTimes] = functionTDT2MatJitterAlign(tdtTimesRaw,mbedTimesRaw);

disp('INITIATING JITTER ALIGNMENT FUNCTION')
%parameters I want to set
bigLim = 2000; %limit for considering something to be excessively long
massDiffNum = 10; %limit to number of big pauses.
medLim = 600; %limit for considering a skipped pulse
expectedITI = 500; %expected ITI

%mbed trace
inputPhotDiff = diff(mbedTimesRaw);

%tdt trace
traceJittDiff = diff(tdtTimesRaw);


%% clean up the MBED signal by removing the big pauses from the checking stage.
disp('Checking for Big Differences in MBED Signal')
findBig = find(inputPhotDiff > bigLim);
%now we need to screen the big differences.
whileCounter = 1;
while length(findBig) >= whileCounter;
    bigSize = findBig(whileCounter);
    if bigSize > length(mbedTimesRaw)/2 %in the case of something coming near the end
        disp('Late Big Diff, deleting')
        mbedTimesRaw(bigSize:end) = [];
        inputPhotDiff = diff(mbedTimesRaw);
        findBig = find(inputPhotDiff > medLim);
    else
        disp('Early Big Difference in Jitter')
    end
    whileCounter = whileCounter + 1;
end

bigDiffs = inputPhotDiff(findBig);

%look for huge diffs. These should indicate the start of the actual session
massDiff = find(bigDiffs > bigLim);
if length(massDiff) <= massDiffNum & length(massDiff)>0
    disp('Long differences in jittered trace. taking last.')
    mbedTimesRaw(1:findBig(massDiff(end))) = [];
    inputPhotDiff = diff(mbedTimesRaw);
elseif isempty(massDiff)
    disp('No long differences in jittered trace')
else
    error('Excessive long TTL Differences in Jittered Trace')
end

%% Now try and fix any skipped TTL pulses that are in the TDT but not the MBED
disp('Looking For Skipped TTL Pulses in MBED Signal')
%find big differences
findBig = find(inputPhotDiff > medLim);
bigDiffs = inputPhotDiff(findBig);
%approximate the length of the differences by expected size. round up. 
bigDiffDiv = round(bigDiffs/expectedITI);
disp('Skipped Pulses Found:')
disp(num2str(length(bigDiffs)))

%now replace with fake points. 
disp('Replacing Failures in MBED Signal')
for i = length(bigDiffs):-1:1
    targetInd = findBig(i)+1;
    mbedTimesRaw(targetInd+(bigDiffDiv-1):end+(bigDiffDiv-1)) = mbedTimesRaw(targetInd:end);
    mbedTimesRaw(targetInd) = mbedTimesRaw(targetInd-1)+expectedITI;
    disp(strcat('Replacing Failure at:',num2str(findBig(i))))
end
inputPhotDiff = diff(mbedTimesRaw);

%% Check that things are aligned with crosscorr
disp('Checking Failure Repair with Cross Correlation')
%now check with crosscorr
[xcf,lags,bounds]  = crosscorr(inputPhotDiff,traceJittDiff,300);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);

if xcLag ~= 0
    disp('CorrLag')
    disp(xcLag)
    disp('MaxCorr')
    disp(xcMax)
    if xcMax < 0.9
        error('Jitter Not Aligned')
    else
        disp('Trying To Fix Jitter With Corr Value')
        if xcLag > 0
            tdtTimesRaw(1:xcLag) = [];
            traceJittDiff = diff(tdtTimesRaw);
            disp('Removed Points from TDT signal')
        elseif xcLag < 0
            mbedTimesRaw(1:-xcLag) = [];
            inputPhotDiff = diff(mbedTimesRaw);
            disp('Removed Points from MBED signal')
        end
    end
elseif xcLag == 0
    disp('Jitter Signal Properly Aligned')
end

%% Confirm that things are okay by looking at lengths.
%now check lengths
if length(inputPhotDiff) ~= length(traceJittDiff)
    disp('Lengths of Jittered traces unequal, attempting removal')
    if length(inputPhotDiff)>length(traceJittDiff) %too many mbed inputs
        mbedTimesRaw(end) = [];
        inputPhotDiff = diff(mbedTimesRaw);
    elseif length(inputPhotDiff)<length(traceJittDiff) %too many TDT inputs
        tdtTimesRaw(length(mbedTimesRaw)+1:end) = [];
        traceJittDiff = diff(tdtTimesRaw);
    end
elseif length(inputPhotDiff) == length(traceJittDiff)
    disp('Lengths of Jittered Traces Equal! YAY')
end

if length(inputPhotDiff) ~= length(traceJittDiff)
    error('Jittered traces still not the right length')
else
    disp('Jittered Traces Now Correct Length')
end

tdtTimes = tdtTimesRaw;
mbedTimes = mbedTimesRaw;



end