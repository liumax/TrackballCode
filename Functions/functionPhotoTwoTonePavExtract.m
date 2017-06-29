%This is a function to extract data regarding basic two tone pavlovian
%task 

function [s] = functionPhotoTwoTonePavExtract(fileName,rasterWindow,thresh,peakToggle);
%% Parameters

% rasterWindow = [-3,5]; %raster window in seconds
% thresh = 0.01; %threshold for peak detection
locoTimeStep = 0.1;

%% Set up structured array for storage
s = struct;

%% First, lets pull the MBED stuff

[trialStates, portStates, trialParams] = maxTrialVariablesLickingTask(fileName);

%we want to look at instates. Port 1 is TDT for photometry, port 2 is
%NOLDUS

inputPhot = [portStates.tStamps',portStates.inStates(:,8)];
outputPhot = [portStates.tStamps',portStates.outStates(:,8)];
rewOut = [portStates.tStamps',portStates.outStates(:,2)];

%eliminate duplicate values!

[inputPhot] = functionSignalDuplicateElim(inputPhot,2);
[outputPhot] = functionSignalDuplicateElim(outputPhot,2);
[rewOut] = functionSignalDuplicateElim(rewOut,2);

%now lets pull some information from the reward outputs. 
rewTimes = rewOut(rewOut(:,2) == 1,1); %pull times
rewDur = rewOut(find(rewOut(:,2) == 1) + 1,1)-rewOut(rewOut(:,2) == 1,1);

%pull just the onset times for photometry
onsetPhot = outputPhot(outputPhot(:,2) == 1,1);
onsetPhotDiff = diff(onsetPhot);

%now lets try and match the trials with rewards
if length(rewTimes) ~= length(onsetPhot)
    toneRewInd = zeros(length(onsetPhot),2);
    toneRewInd(:,1) = 1:length(onsetPhot);
    for i = 1:length(onsetPhot)
        testDiff = abs(rewTimes - onsetPhot(i));
        if min(testDiff) < 5000
            [minVal minInd] = min(testDiff);
    %         toneRewInd(i,1) = i;
            toneRewInd(i,2) = minInd;
        end
    end

    trialLow = find(toneRewInd(:,2) == 0);
    trialHi = find(toneRewInd(:,2) ~= 0);
else
    trialLow = find(rewDur == min(rewDur));
    trialHi = find(rewDur == max(rewDur));
    toneRewInd = zeros(length(onsetPhot),2);
    toneRewInd(:,1) = 1:length(onsetPhot);
    toneRewInd(:,2) = 1:length(onsetPhot);
end


%Now pull locomotor data

[locoData] = functionMBEDrotary(portStates.inStates(:,4),portStates.inStates(:,5),portStates.tStamps/1000,locoTimeStep);

%% Now lets pull the photometry inputs

%load file
data = load(fileName);
data=data.data;


%Code from Chris that performs isosbestic correction. 
traceDF = isosbestic_correction(data);

%pull timestamps for fluorescence
traceTiming = [0:1/data.streams.x70G.fs:(1/data.streams.x70G.fs)*(length(data.streams.x70G.data)-1)];


%pull peaks 170616 This appears to have problem: built for 2016 matlab, has
%additional functionality for peak finding.
if peakToggle == 1
    [peakInfo, riseInfo, troughInfo] = findPhotoPeaks(traceTiming,traceDF,thresh);
end

%pull jittered signal
traceJitt = data.epocs.PtE1.onset;
traceJittDiff = diff(traceJitt);
%pull output timings
traceMBED = data.epocs.PtC0.onset;
traceMBEDDiff = diff(traceMBED);

%check alignment
[xcf,lags,bounds]  = crosscorr(onsetPhotDiff,traceMBEDDiff);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);

if xcLag ~= 0
    error('Photometry Not Aligned')
elseif length(onsetPhot) ~= length(traceMBED)
    error('Mismatch in Number of Photometry Pulses')
end

%calculate raster in terms of time steps in photometry
photoTimeStep = 1/data.streams.x70G.fs;
rasterPhotWindow = round(rasterWindow/photoTimeStep);

%now lets get the signal inputs cleaned up! This should eliminate sets of
%extra pulses from both before and after the actual recording session.
inputPhotOnset = inputPhot(inputPhot(:,2) == 1,1);
inputPhotOnsetDiff = diff(inputPhotOnset);

bigDiffFind = find(inputPhotOnsetDiff>550);
while length(bigDiffFind)>0
    targetInd = bigDiffFind(1);
    if targetInd < length(traceJitt)
%         inputPhot = inputPhot(targetInd*2+1:end,:);
        inputPhotOnset = inputPhotOnset(targetInd+1:end);
        inputPhotOnsetDiff = diff(inputPhotOnset);
        bigDiffFind = find(inputPhotOnsetDiff>550);
    elseif targetInd >= length(traceJitt)
%         inputPhot = inputPhot(1:targetInd*2+1);
        inputPhotOnset = inputPhotOnset(1:targetInd);
        inputPhotOnsetDiff = diff(inputPhotOnset);
        bigDiffFind = find(inputPhotOnsetDiff>550);
    end
end

%okaaaay I think I have an excess datapoint possibly. Lets build out code
%to eliminate this. First, do crosscorr to get proper alignment. Then trim
%off the end.

[xcf,lags,bounds] = crosscorr(traceJittDiff,inputPhotOnsetDiff);
[MAX maxInd] = max(xcf);
maxLag = lags(maxInd);


while maxLag ~= 0
    disp(maxLag)
    disp('Jittered Traces Not Aligned, Adjusting')
    if maxLag > 0
        inputPhotOnset = inputPhotOnset(1+maxLag:end);
        inputPhotOnsetDiff = diff(inputPhotOnset);
    elseif maxLag < 0
        traceJitt = traceJitt(1-maxLag:end);
        traceJittDiff = diff(traceJitt);
    end
    [xcf,lags,bounds] = crosscorr(traceJittDiff,inputPhotOnsetDiff);
    [MAX maxInd] = max(xcf);
    maxLag = lags(maxInd);
end

disp('Traces Aligned')

while length(traceJittDiff) ~= length(inputPhotOnsetDiff)
    if length(traceJittDiff) > length(inputPhotOnsetDiff)
        diffFind = length(traceJittDiff) - length(inputPhotOnsetDiff);
        traceJitt = traceJitt(1:end-diffFind);
        traceJittDiff = diff(traceJitt);
        disp('Lengths Adjusted by Shortening TraceJitt')
    elseif length(traceJittDiff) < length(inputPhotOnsetDiff)
        diffFind = length(traceJittDiff) - length(inputPhotOnsetDiff);
        inputPhotOnset = inputPhotOnset(1:end+diffFind);
        inputPhotOnsetDiff = diff(inputPhotOnset);
        disp('Lengths Adjusted by Shortening InputPhot')
    end
end



%% Now lets sort out responses

%First, check to see if number of tones equals number of outputs.

if length(onsetPhot) == length(traceMBED)
    disp('Trial Number and Photometry Inputs Matched')
elseif length(onsetPhot) ~= length(traceMBED)
    error('Mismatch in Trial Number and Photometry Inputs')
end

%then, lets just pull the basic raster

photoRaster = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,length(traceMBED));

for i = 1:length(traceMBED)
    alignTime = traceMBED(i);
    %find the time in the photometry trace
    photoPoint = find(traceTiming - alignTime > 0,1,'first');
    if photoPoint + rasterPhotWindow(2) < length(traceDF)
        photoRaster(:,i) = traceDF(photoPoint + rasterPhotWindow(1):photoPoint + rasterPhotWindow(2));
    else
        break
    end
end


photoRasterRew = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,length(rewTimes));

for i = 1:length(rewTimes)
    alignTime = traceMBED(i);
    %find the time in the photometry trace
    photoPoint = find(traceTiming - alignTime > 0,1,'first');
    photoRasterRew(:,i) = traceDF(photoPoint + rasterPhotWindow(1):photoPoint + rasterPhotWindow(2));
end

%store means of different tones
meanHi = mean(photoRaster(:,trialHi)');
meanLow = mean(photoRaster(:,trialLow)');
steHi = std(photoRaster(:,trialHi)')/sqrt(length(trialHi));
steLow = std(photoRaster(:,trialLow)')/sqrt(length(trialHi));


%% Velocity Data
% Here, we will use the photometry as the base, since it is probably more
% reliable

velTrueTime = interp1(onsetPhot/1000,traceMBED,locoData.Velocity(:,1));

%in the current iteration, this has issues because the photometry cant
%align without pulses, and i only have output pulses from the MBED when the
%sound turns on. Therefore, I will need to crop things out.

findVelFirst = find(~isnan(velTrueTime),1,'first');
findVelLast = find(~isnan(velTrueTime),1,'last');
%find nearest in photometry signal
findPhotFirst = find(traceTiming - velTrueTime(findVelFirst)>0,1,'first');
findPhotLast = find(traceTiming - velTrueTime(findVelLast)>0,1,'first');

%now lets align the velocity with tone
%since step size for locomotion is set, lets use that to calculate raster
%windows
velWindow = rasterWindow/locoTimeStep;
velRasterAxis = [rasterWindow(1):locoTimeStep:rasterWindow(2)];
for i = 1:length(traceMBED)
    targetTime = traceMBED(i);
    findTime = find(velTrueTime - targetTime > 0,1,'first');
    if findTime + velWindow(2) < length(velTrueTime) & findTime + velWindow(1) > 0
        velRaster(:,i) = locoData.Velocity((findTime + velWindow(1)):(findTime + velWindow(2)),2);
    else
        disp('EDGE ISSUE')
        disp(i)
        velRaster(:,i) = zeros(length(velRasterAxis),1);
    end
end

%get means

velHiMean = mean(velRaster(:,trialHi)');
velHiSTE = std(velRaster(:,trialHi)')/sqrt(length(trialHi));

velLowMean = mean(velRaster(:,trialLow)');
velLowSTE = std(velRaster(:,trialLow)')/sqrt(length(trialLow));

%now do the same for reward times
velWindowRew = rasterWindow/locoTimeStep;
velRasterAxisRew = [rasterWindow(1):locoTimeStep:rasterWindow(2)];
for i = 1:length(rewTimes)
    targetTime = rewTimes(i)/1000;
    findTime = find(locoData.Velocity(:,1) - targetTime > 0,1,'first');
    if findTime + velWindow(2) < length(locoData.Velocity(:,1)) & findTime + velWindow(1) > 0
        velRasterRew(:,i) = locoData.Velocity((findTime + velWindow(1)):(findTime + velWindow(2)),2);
    else
        disp('EDGE ISSUE')
        disp(i)
        velRasterRew(:,i) = zeros(length(velRasterAxis),1);
    end
end


%% Licking Data

lickData = trialParams.licking;
if length(lickData) > 0
    lickTrueTimes = interp1(onsetPhot/1000,traceMBED,lickData(:,1)/1000);
    
    [lickRasterRew] = functionBasicRaster(lickData(:,1)/1000,rewTimes/1000,rasterWindow);
    [lickRasterTone] = functionBasicRaster(lickData(:,1)/1000,onsetPhot/1000,rasterWindow);
    %now match to high/low
    lickInd = 1;
    for i = 1:length(trialHi)
        %see if there are values that are present
        lickFinder = find(lickRasterTone(:,2) == trialHi(i));
        if length(lickFinder) > 0
            lickRasterToneHi(lickInd:lickInd + length(lickFinder)-1,:) = lickRasterTone(lickFinder,:);
            lickInd = lickInd + length(lickFinder);
        end
    end

    lickInd = 1;
    for i = 1:length(trialLow)
        %see if there are values that are present
        lickFinder = find(lickRasterTone(:,2) == trialLow(i));
        if length(lickFinder) > 0
            lickRasterToneLow(lickInd:lickInd + length(lickFinder)-1,:) = lickRasterTone(lickFinder,:);
            lickInd = lickInd + length(lickFinder);
        end
    end

    %store as histograms
    histLickToneHi = hist(lickRasterToneHi(:,1),[rasterWindow(1):0.1:rasterWindow(2)]);
    histLickToneLow = hist(lickRasterToneLow(:,1),[rasterWindow(1):0.1:rasterWindow(2)]);
    histLickRewHi = hist(lickRasterRew(:,1),[rasterWindow(1):0.1:rasterWindow(2)]);
    histLickAxis = [rasterWindow(1):0.1:rasterWindow(2)];
else
    lickTrueTimes = [];
    lickRasterRew = [0,0];
    lickRasterTone = [0,0];
    histLickAxis = [rasterWindow(1):0.1:rasterWindow(2)];
    histLickToneHi = zeros(length(histLickAxis),1);
    histLickToneLow = zeros(length(histLickAxis),1);
    histLickRewHi = zeros(length(histLickAxis),1);
end



%% Store things into the structure
s.Trials = toneRewInd;

s.Velocity.RasterAxis = velRasterAxis;
s.Velocity.ToneRaster = velRaster;
s.Velocity.RewardRaster = velRasterRew;

s.Licking.ToneRaster = lickRasterTone;
s.Licking.RewardRaster = lickRasterRew;


%create a set of points for the raster axis
rasterAxis = zeros(rasterWindow(2)-rasterWindow(1)+1,2);
rasterAxis(:,1) = [rasterWindow(1):rasterWindow(2)];
%find a second in photometry sample time
rasterAxis(:,2) = [0:length(photoRaster)/(rasterWindow(2) - rasterWindow(1)):length(photoRaster)];
rasterAxis(1,2) = 1;

s.Photometry.RasterAxis = rasterAxis;
s.Photometry.ToneRaster = photoRaster;
s.Photometry.RewardRaster = photoRasterRew;

end





















