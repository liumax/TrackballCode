%This is code that will read and link MBED inputs with photometry, while
%recording from my big rig. 
function [s] = analysisPavTwoToneCatchFree(fileName);
%% Parameters

rasterWindow = [-3,5]; %raster window in seconds
thresh = 0.01; %threshold for peak detection
locoTimeStep = 0.1;
minBar = 8000; %for clearing out false pulses from speaker card.

%store things in a big structure
s = struct;
s.Parameters.RasterWindow = rasterWindow;
s.Parameters.PeakThreshold = thresh;
s.Parameters.LocomotionTimeStep = locoTimeStep;

lickHistBin = 0.1;
%% First, lets pull the MBED stuff


[trialStates, portStates, trialParams] = maxTrialVariablesLickingTask(fileName);

%we want to look at instates. Port 1 is TDT for photometry, port 2 is
%NOLDUS

inputPhot = [portStates.tStamps',portStates.inStates(:,8)];
outputPhot = [portStates.tStamps',portStates.outStates(:,8)];
rewOut = [portStates.tStamps',portStates.outStates(:,2)];
toneIn = [portStates.tStamps',portStates.inStates(:,1)];

%eliminate duplicate values!
[inputPhot] = functionSignalDuplicateElim(inputPhot,2);
[outputPhot] = functionSignalDuplicateElim(outputPhot,2);
[rewOut] = functionSignalDuplicateElim(rewOut,2);
[toneIn] = functionSignalDuplicateElim(toneIn,2);

%check the toneIn, separate out to both onset times and duration
toneOnset = toneIn(toneIn(:,2) == 1,1);
toneFirst = find(toneIn(:,2) == 1,1,'first');
toneDurs = toneIn([toneFirst+1:2:end],1) - toneIn([toneFirst:2:end],1);

%now lets pull some information from the reward outputs. 
rewTimes = rewOut(rewOut(:,2) == 1,1); %pull times
rewDur = rewOut(find(rewOut(:,2) == 1) + 1,1)-rewOut(rewOut(:,2) == 1,1);

%pull just the onset times for photometry
onsetPhot = outputPhot(outputPhot(:,2) == 1,1);
firstFind = find(outputPhot(:,2) == 1,1,'first');
outputDur = (outputPhot(firstFind+1:2:end,1))-(outputPhot(firstFind:2:end,1));
onsetPhotDiff = diff(onsetPhot);

inputPhotOnset = inputPhot(inputPhot(:,2) == 1,1);
inputPhotDiff = diff(inputPhotOnset);

%toggle for marking when there are MBED issues with outputs. 
mbedErrorFlag = 0;

%171011, due to experience with this and some reliability issues I've
%observed, will use the input (toneOnset) rather than the output
%(outputPhot) for determining tone timing. 

%check for errors by comparing trialStates with actual tone inputs
if length(trialStates.trialType) ~= length(toneOnset)
    error('INCORRECT NUMBER OF TRIALS DETECTED. ABORT')
end

%now check for the different trial types available
disp(strcat('Low Tone Trials:',num2str(length(find(trialStates.trialType == 1)))))
disp(strcat('Hi Tone Trials:',num2str(length(find(trialStates.trialType == 2)))))
disp(strcat('Free Rew Trials:',num2str(length(find(trialStates.trialType == 3)))))
disp(strcat('Catch Tone Trials:',num2str(length(find(trialStates.trialType == 4)))))

%now save these trial types as variables
trialLow = find(trialStates.trialType == 1);
trialHi = find(trialStates.trialType == 2);
trialFree = find(trialStates.trialType == 3);
trialCatch = find(trialStates.trialType == 4);
trialRew = sort([trialHi,trialFree]);

%NOW DETERMINE the difference between the outputs (outputPhot) and the
%actual inputs
if length(onsetPhot) == length(toneOnset)
    disp('Tone Input and Output Numbers Matched, Proceed!')
    delInds = [];
else
    disp('Tone Input and Output Numbers Mismatched, Correcting...')
    %calculate difference
    pulseDiff = length(onsetPhot) - length(toneOnset);
    disp(strcat('Output has ',num2str(pulseDiff),' extra pulses'))
    %find intersect of the two
    [C,ia,ib] = intersect(onsetPhot,toneOnset); %ia is the index for onsetPhot. 
    delInds = [1:length(onsetPhot)];
    delInds(ia) = [];
    %check to make sure this is correct!
    if length(delInds) ~= pulseDiff
        error('FAILURE TO CORRECT OUTPUT PROPERLY')
    end
    onsetPhot(delInds) = [];
    
    %confirm that things line up properly!
    [xcf,lags,bounds]  = crosscorr(onsetPhot,toneOnset);
    [xcMax maxInd] = max(xcf);
    xcLag = lags(maxInd);
    
    if xcLag ~= 0
        error('Tone Inputs Not Aligned to Tone Outputs!')
    end
    
end


s.MBED.RewTimes = rewTimes;
s.MBED.RewDur = rewDur;
% s.MBED.Jitter = inputPhotOnset;
s.MBED.ToneDelivery = toneOnset;
s.MBED.ToneOutput = onsetPhot;
s.MBED.HiTrials = trialHi;
s.MBED.LowTrials = trialLow;
s.MBED.FreeTrials = trialFree;
s.MBED.CatchTrials = trialCatch;
s.MBED.RewTrials = trialRew;
s.MBED.Raw = portStates;
s.MBED.ToneOutputFailures = delInds;
%Now pull locomotor data

[locoData] = functionMBEDrotary(portStates.inStates(:,4),portStates.inStates(:,5),portStates.tStamps/1000,locoTimeStep);

s.Locomotion = locoData;

%% Now lets pull the photometry inputs

%load file
data = load(strcat(fileName,'.mat'));
data=data.data;


% %Code from Chris that performs isosbestic correction. 
% traceDF = isosbestic_correction(data);
% 
% 
% %pull timestamps for fluorescence
% traceTiming = [0:1/data.streams.x70G.fs:(1/data.streams.x70G.fs)*(length(data.streams.x70G.data)-1)];
% 
% s.Photo.dFTrace = traceDF;
% s.Photo.dFTime = traceTiming;
% s.Photo.x70 = data.streams.x70G.data;
% s.Photo.x05 = data.streams.x05G.data;
% s.Photo.Raw = data.streams.Fi1r.data;
% s.Photo.RawRate = data.streams.Fi1r.fs;


%170809 Replacing with scott code
[filtSig1,filtSig2,traceDF,traceTiming] = functionPhotometryRawExtraction(data);

s.Photo.dFTrace = traceDF;
s.Photo.dFTime = traceTiming;
s.Photo.x70 = filtSig1;
s.Photo.x05 = filtSig2;
s.Photo.Raw = data.streams.Fi1r.data;
s.Photo.RawRate = data.streams.Fi1r.fs;


%pull peaks 170616 This appears to have problem: built for 2016 matlab, has
%additional functionality for peak finding.
try
    [t_ds,newSmoothDS,targetPeaks] = functionPhotoPeakProcess(traceTiming,filtSig1,0.01);
%     [peakInfo, riseInfo, troughInfo] = findPhotoPeaks(traceTiming,traceDF,thresh);
catch
    disp('Peak Detection Failed')
    targetPeaks = [];
    newSmoothDS = [];
    t_ds = [];
end

s.Photo.Peaks = targetPeaks;
s.Photo.Photo.x70dF = newSmoothDS;
s.Photo.Photo.x70dFTime = t_ds;

%pull jittered signal
traceJitt = data.epocs.PtE1.onset;
traceJittDiff = diff(traceJitt);

%clean up the MBED signal. generally there will be excess TTLs at the
%beginning
findBig = find(inputPhotDiff > 2000);
%now we need to screen the big differences.
whileCounter = 1;
while length(findBig) >= whileCounter;
    bigSize = findBig(whileCounter);
    if bigSize > length(inputPhotOnset)/2 %in the case of something coming near the end
        disp('Late Big Diff, deleting')
        inputPhotOnset(bigSize:end) = [];
        inputPhotDiff = diff(inputPhotOnset);
        findBig = find(inputPhotDiff > 600);
    else
        disp('Early Big Difference in Jitter')
    end
    whileCounter = whileCounter + 1;
end

bigDiffs = inputPhotDiff(findBig);

%look for huge diffs. These should indicate the start of the actual session
massDiff = find(bigDiffs > 2000);
if length(massDiff) <=3 & length(massDiff)>0
    disp('Long Differences in jittered trace. taking last.')
    inputPhotOnset(1:findBig(massDiff(end))) = [];
    inputPhotDiff = diff(inputPhotOnset);
elseif isempty(massDiff)
    disp('No Long differences in jittered trace')
else
    error('Excessive Large TTL Differences in Jittered Trace')
end

%find big differences
findBig = find(inputPhotDiff > 600);
bigDiffs = inputPhotDiff(findBig);
%approximate the length of the differences by 500. round up. 
bigDiffDiv = round(bigDiffs/500);

%now replace with fake points. 
for i = length(bigDiffs):-1:1
    targetInd = findBig(i)+1;
    inputPhotOnset(targetInd+(bigDiffDiv-1):end+(bigDiffDiv-1)) = inputPhotOnset(targetInd:end);
    inputPhotOnset(targetInd) = inputPhotOnset(targetInd-1)+500;
end
inputPhotDiff = diff(inputPhotOnset);

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
            traceJitt(1:xcLag) = [];
            traceJittDiff = diff(traceJitt);
        elseif xcLag < 0
            inputPhotOnset(1:-xcLag) = [];
            inputPhotDiff = diff(inputPhotOnset);
        end
    end
elseif xcLag == 0
    disp('Jitter Signal Properly Aligned')
end



%now check lengths
if length(inputPhotDiff) ~= length(traceJittDiff)
    disp('Lengths of Jittered traces unequal, attempting removal')
    if length(inputPhotDiff)>length(traceJittDiff) %too many mbed inputs
        inputPhotOnset(end) = [];
        inputPhotDiff = diff(inputPhotOnset);
    elseif length(inputPhotDiff)<length(traceJittDiff)
        traceJitt(length(inputPhotOnset)+1:end) = [];
        traceJittDiff = diff(traceJitt);
    end
elseif length(inputPhotDiff) == length(traceJittDiff)
    disp('Lengths of Jittered Traces Equal! YAY')
end

if length(inputPhotDiff) ~= length(traceJittDiff)
    error('Jittered traces still not the right length')
else
    disp('Jittered Traces Now Correct Length')
end

s.Photo.Jitter = traceJitt;
s.MBED.Jitter = inputPhotOnset;

%now lets deal with inputs. This is now designed to try and handle cases in
%which sound pulses are not delivered. 
% pull output timings
try
    traceMBED = data.epocs.PtC0.onset;
    interpTrig = 0;
catch
    disp('No TDT Tone Pulses Detected: Using MBED Tone Times to Interpolate')
    traceMBED = interp1(inputPhotOnset,traceJitt,toneOnset);
    interpTrig = 1;
end
traceMBEDDiff = diff(traceMBED);

%stash a copy of deletion indices.
delBackup = delInds;

%check alignment
if mbedErrorFlag == 1
    disp('Correcting False Signals in TDT MBED SIGNAL')
    traceMBED(delInds) = [];
end
traceMBEDDiff = diff(traceMBED);

[xcf,lags,bounds]  = crosscorr(onsetPhotDiff/1000,traceMBEDDiff);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);



if xcLag ~= 0
    disp('CorrLag')
    disp(xcLag)
    disp('MaxCorr')
    disp(xcMax)
    if xcMax < 0.9
        error('Tones Not Aligned')
    else
        disp('Trying To Fix Tone Signal With Corr Value')
        if xcLag > 0
            traceMBED(1:xcLag) = [];
            traceMBEDDiff = diff(traceMBED);
        elseif xcLag < 0
            onsetPhot(1:-xcLag) = [];
            onsetPhotDiff = diff(onsetPhot);
        end
    end
elseif xcLag == 0
    disp('Tone Signal Properly Aligned')
end

%now check lengths
if length(onsetPhotDiff) ~= length(traceMBEDDiff)
    disp(strcat('OnsetPhot',num2str(length(onsetPhot)),'TraceMBED:',num2str(length(traceMBED))))
    error('Tone Signal Hasnt Been Corrected')
elseif length(onsetPhotDiff) == length(traceMBEDDiff)
    disp('Lengths of Tone Signals Equal! YAY')
end

s.Photo.MBEDSig = traceMBED; %store tone times, makes life easier later on. 

%calculate raster in terms of time steps in photometry
photoTimeStep = mean(diff(t_ds));
rasterPhotWindow = round(rasterWindow/photoTimeStep);
rasterVect = [rasterWindow(1):photoTimeStep:rasterWindow(2)];

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
    photoPoint = find(t_ds - alignTime > 0,1,'first');
    if photoPoint + rasterPhotWindow(2) < length(newSmoothDS)
        photoRaster(:,i) = newSmoothDS(photoPoint + rasterPhotWindow(1):photoPoint + rasterPhotWindow(2));
    else
        disp('Tone Rasters: Reached End of Photometry Trace')
        disp(i)
        break
    end
end

%now lets try and find things!
zeroPoint = find(rasterVect > 0,1,'first');
preBin = -10;
postBin = 20;
baseVal = min(photoRaster(zeroPoint + preBin:zeroPoint,:));
peakVal = max(photoRaster(zeroPoint:zeroPoint + postBin,:));
magVal = peakVal - baseVal;

valStore = [baseVal;peakVal;magVal];

s.Vals = valStore;
%align to reward times
if length(rewTimes)>0
    photoRasterRew = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,length(rewTimes));

    %generate the correct times by interpolation
    rewTimeTDT = interp1(inputPhotOnset,traceJitt,rewTimes);

    for i = 1:length(rewTimes)
        alignTime = rewTimeTDT(i);
        %find the time in the photometry trace
        photoPoint = find(t_ds - alignTime > 0,1,'first');
        if photoPoint + rasterPhotWindow(2) < length(newSmoothDS)
            photoRasterRew(:,i) = newSmoothDS(photoPoint + rasterPhotWindow(1):photoPoint + rasterPhotWindow(2));
        else
            disp('Rew Rasters: Reached End of Photometry Trace')
            disp(i)
            break
        end
    end
else
    photoRasterRew = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,1)
end

%now try aligning to licks
if length(trialParams.licking)>0
    lickTimesTDT = interp1(inputPhotOnset,traceJitt,trialParams.licking(:,1));
    %lets use smaller raster for licking
    lickRastWindow = [-1 1];
    lickPhotWindow = round(lickRastWindow/photoTimeStep);

    for i = 1:length(rewTimes)
        alignTime = rewTimeTDT(i);
        %find the time in the photometry trace
        photoPoint = find(t_ds - alignTime > 0,1,'first');
        if photoPoint + rasterPhotWindow(2) < length(newSmoothDS)
            photoRasterLick(:,i) = newSmoothDS(photoPoint + lickPhotWindow(1):photoPoint + lickPhotWindow(2));
        else
            disp('LickRasters: Reached End of Photometry Trace')
            disp(i)
            break
        end
    end
else
    %lets use smaller raster for licking
    lickRastWindow = [-1 1];
    lickPhotWindow = round(lickRastWindow/photoTimeStep);
    photoRasterLick = zeros(lickPhotWindow(2) - lickPhotWindow(1)+1,1);
end



%store means of different tones



steHi = std(photoRaster(:,trialHi)')/sqrt(length(trialHi));
steLow = std(photoRaster(:,trialLow)')/sqrt(length(trialLow));
steFree = std(photoRaster(:,trialFree)')/sqrt(length(trialFree));
steCatch = std(photoRaster(:,trialCatch)')/sqrt(length(trialCatch));
steRew = std(photoRasterRew')/sqrt(length(traceMBED));

meanHi = [mean(photoRaster(:,trialHi)');mean(photoRaster(:,trialHi)') - steHi;mean(photoRaster(:,trialHi)') + steHi];
meanLow = [mean(photoRaster(:,trialLow)');mean(photoRaster(:,trialLow)') - steLow;mean(photoRaster(:,trialLow)') + steLow];
meanFree = [mean(photoRaster(:,trialFree)');mean(photoRaster(:,trialFree)') - steLow;mean(photoRaster(:,trialFree)') + steFree];
meanCatch = [mean(photoRaster(:,trialCatch)');mean(photoRaster(:,trialCatch)') - steCatch;mean(photoRaster(:,trialCatch)') + steCatch];
meanRew = [mean(photoRasterRew');mean(photoRasterRew') - steRew;mean(photoRasterRew') + steRew];

photoZero = (-1*rasterWindow(1))/photoTimeStep;

s.PhotoRaster.ToneRaster = photoRaster;
s.PhotoRaster.RewardRaster = photoRasterRew;
s.PhotoRaster.MeanRew = meanRew;
s.PhotoRaster.MeanHi = meanHi;
s.PhotoRaster.MeanLow = meanLow;
s.PhotoRaster.MeanFree = meanFree;
s.PhotoRaster.MeanCatch = meanCatch;
s.PhotoRaster.LickRaster = photoRasterLick;
s.PhotoRaster.TimeStep = photoTimeStep;


%% Velocity Data
% Here, we will use the photometry as the base, since it is probably more
% reliable
if length(locoData.Velocity) > 0
    velTrueTime = interp1(onsetPhot/1000,traceMBED,locoData.Velocity(:,1));

    %in the current iteration, this has issues because the photometry cant
    %align without pulses, and i only have output pulses from the MBED when the
    %sound turns on. Therefore, I will need to crop things out.

    findVelFirst = find(~isnan(velTrueTime),1,'first');
    findVelLast = find(~isnan(velTrueTime),1,'last');
    %find nearest in photometry signal
    findPhotFirst = find(traceTiming - velTrueTime(findVelFirst)>0,1,'first');
    findPhotLast = find(traceTiming - velTrueTime(findVelLast)>0,1,'first');
else
    velTrueTime = [];
    findVelFirst = [];
    findVelLast = [];
    findPhotFirst = 1;
    findPhotLast = length(traceTiming);
end

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
        disp('Velocity EDGE ISSUE: tone alignment')
        i
        disp(num2str(length(velTrueTime)-(findTime + velWindow(2))))
%         disp(i)
        velRaster(:,i) = zeros(velWindow(2)-velWindow(1)+1,1);
%         break
    end
end

%get means

velHiMean = mean(velRaster(:,trialHi)');
velHiSTE = std(velRaster(:,trialHi)')/sqrt(length(trialHi));

velLowMean = mean(velRaster(:,trialLow)');
velLowSTE = std(velRaster(:,trialLow)')/sqrt(length(trialLow));

velFreeMean = mean(velRaster(:,trialFree)');
velFreeSTE = std(velRaster(:,trialFree)')/sqrt(length(trialFree));

velCatchMean = mean(velRaster(:,trialCatch)');
velCatchSTE = std(velRaster(:,trialCatch)')/sqrt(length(trialCatch));

%now do the same for reward times
velWindowRew = rasterWindow/locoTimeStep;
velRasterAxisRew = [rasterWindow(1):locoTimeStep:rasterWindow(2)];
for i = 1:length(rewTimes)
    targetTime = rewTimes(i)/1000;
    findTime = find(locoData.Velocity(:,1) - targetTime > 0,1,'first');
    if findTime + velWindow(2) < length(locoData.Velocity(:,1)) & findTime + velWindow(1) > 0
        velRasterRew(:,i) = locoData.Velocity((findTime + velWindow(1)):(findTime + velWindow(2)),2);
    else
        disp('Velocity EDGE ISSUE, Reward Alignment')
%         disp(i)
        velRasterRew(:,i) = zeros(velWindow(2)-velWindow(1)+1,1);
%         break
    end
end

velRewMean = mean(velRasterRew');
velRewSTE = std(velRasterRew')/sqrt(length(trialRew));

s.VelRaster.ToneRaster = velRaster;
s.VelRaster.RewardRaster = velRasterRew;
s.VelRaster.Axis = velRasterAxis;
s.VelRaster.MeanRew = [velRewMean;velRewMean-velRewSTE;velRewMean+velRewSTE];
s.VelRaster.MeanHi = [velHiMean;velHiMean-velHiSTE;velHiMean+velHiSTE];
s.VelRaster.MeanLow = [velLowMean;velLowMean-velLowSTE;velLowMean+velLowSTE];
s.VelRaster.MeanFree = [velFreeMean;velFreeMean-velFreeSTE;velFreeMean+velFreeSTE];
s.VelRaster.MeanCatch = [velCatchMean;velCatchMean-velCatchSTE;velCatchMean+velCatchSTE];

%perform ROC analysis
% [funcOut] = functionROCLocoPhoto(s.MBED.Jitter,s.Photo.Jitter,s.Photo.Photo.x70dF,s.Photo.Photo.x70dFTime,s.Locomotion.Velocity);
% s.ROC = funcOut;
%% Licking Data

lickData = trialParams.licking;
s.MBED.Licks = lickData;

if length(lickData)>0

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
    lickRasterToneLow = [0 0];
    for i = 1:length(trialLow)
        %see if there are values that are present
        lickFinder = find(lickRasterTone(:,2) == trialLow(i));
        
        if length(lickFinder) > 0
            lickRasterToneLow(lickInd:lickInd + length(lickFinder)-1,:) = lickRasterTone(lickFinder,:);
            lickInd = lickInd + length(lickFinder);
        end
    end
    
    lickInd = 1;
    lickRasterToneFree = [0 0];
    for i = 1:length(trialFree)
        %see if there are values that are present
        lickFinder = find(lickRasterTone(:,2) == trialFree(i));
        
        if length(lickFinder) > 0
            lickRasterToneFree(lickInd:lickInd + length(lickFinder)-1,:) = lickRasterTone(lickFinder,:);
            lickInd = lickInd + length(lickFinder);
        end
    end
    
    lickInd = 1;
    lickRasterToneCatch = [0 0];
    for i = 1:length(trialCatch)
        %see if there are values that are present
        lickFinder = find(lickRasterTone(:,2) == trialCatch(i));
        
        if length(lickFinder) > 0
            lickRasterToneCatch(lickInd:lickInd + length(lickFinder)-1,:) = lickRasterTone(lickFinder,:);
            lickInd = lickInd + length(lickFinder);
        end
    end

    %store as histograms
    histLickToneHi = hist(lickRasterToneHi(:,1),[rasterWindow(1):lickHistBin:rasterWindow(2)]);
    histLickToneLow = hist(lickRasterToneLow(:,1),[rasterWindow(1):lickHistBin:rasterWindow(2)]);
    histLickToneFree = hist(lickRasterToneFree(:,1),[rasterWindow(1):lickHistBin:rasterWindow(2)]);
    histLickToneCatch = hist(lickRasterToneCatch(:,1),[rasterWindow(1):lickHistBin:rasterWindow(2)]);
    
    histLickRewHi = hist(lickRasterRew(:,1),[rasterWindow(1):lickHistBin:rasterWindow(2)]);
    histLickAxis = [rasterWindow(1):lickHistBin:rasterWindow(2)];
else
    lickRasterTone = zeros(1,2);
    lickRasterRew = zeros(1,2);
    lickRasterToneHi = zeros(1,2);
    lickRasterToneLow = zeros(1,2);
    histLickToneHi = zeros(length([rasterWindow(1):lickHistBin:rasterWindow(2)]),1);
    histLickToneLow = zeros(length([rasterWindow(1):lickHistBin:rasterWindow(2)]),1);
    histLickToneFree = zeros(length([rasterWindow(1):lickHistBin:rasterWindow(2)]),1);
    histLickToneCatch = zeros(length([rasterWindow(1):lickHistBin:rasterWindow(2)]),1);
    histLickRewHi = zeros(length([rasterWindow(1):lickHistBin:rasterWindow(2)]),1);
    histLickAxis = [rasterWindow(1):lickHistBin:rasterWindow(2)];
end

%bin pre tones. 
%first, find minimum tone reward latency
try
    toneRewLag = min(rewTimes - (onsetPhot(trialHi)));
catch
    disp('Failure to Find ToneRewLag')
    toneRewLag = 1300;
    length(rewTimes)
    length(trialHi)
end
%now find tones by trial
binLick = zeros(length(onsetPhot),1);
for i = 1:length(onsetPhot)
    %find if there is any licks
    lickFinder = find(lickRasterTone(:,2) == i & lickRasterTone(:,1) > 0 & lickRasterTone(:,1) < toneRewLag/1000);
    if lickFinder
        binLick(i) = length(lickFinder);
    end
end


s.Licking.ToneRaster = lickRasterTone/length(toneOnset)/lickHistBin;
s.Licking.RewardRaster = lickRasterRew/length(trialRew)/lickHistBin;
s.Licking.ToneHistHi = histLickToneHi/length(trialHi)/lickHistBin;
s.Licking.ToneHistLow = histLickToneLow/length(trialLow)/lickHistBin;
s.Licking.ToneHistRew = histLickRewHi/length(trialRew)/lickHistBin;
s.Licking.ToneHistFree = histLickToneFree/length(trialFree)/lickHistBin;
s.Licking.ToneHistCatch = histLickToneCatch/length(trialCatch)/lickHistBin;
s.Licking.Axis = histLickAxis;
s.LickingToneBin = binLick;

        
%% Plot everything
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%create a set of points for the raster axis
rasterAxis = zeros(rasterWindow(2)-rasterWindow(1)+1,2);
rasterAxis(:,1) = [rasterWindow(1):rasterWindow(2)];
%find a second in photometry sample time
rasterAxis(:,2) = [0:length(photoRaster)/(rasterWindow(2) - rasterWindow(1)):length(photoRaster)];
rasterAxis(1,2) = 1;

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])

%plot overall photometry trace and locomotion trace
subplot(4,3,1)
hold on
if length(locoData.Velocity > 0)
    plot(velTrueTime(findVelFirst:findVelLast),(locoData.Velocity(findVelFirst:findVelLast,2)-min(locoData.Velocity(findVelFirst:findVelLast,2)))/(max(locoData.Velocity(findVelFirst:findVelLast,2))-min(locoData.Velocity(findVelFirst:findVelLast,2))))
    plot(traceTiming(findPhotFirst:1000:findPhotLast),(traceDF(findPhotFirst:1000:findPhotLast)-min(traceDF(findPhotFirst:findPhotLast)))/(max(traceDF(findPhotFirst:findPhotLast))-min(traceDF(findPhotFirst:findPhotLast))),'r')
    xlim([velTrueTime(findVelFirst),velTrueTime(findVelLast)])
else
    plot(traceTiming(findPhotFirst:1000:findPhotLast),(traceDF(findPhotFirst:1000:findPhotLast)-min(traceDF(findPhotFirst:findPhotLast)))/(max(traceDF(findPhotFirst:findPhotLast))-min(traceDF(findPhotFirst:findPhotLast))),'r')
    xlim([traceTiming(findPhotFirst),traceTiming(findPhotLast)])
end
title('Normalized Vel (b) and Photometry (r)')

%plot average hi and low traces for photometry
subplot(4,3,4)
hold on
plot(meanHi(1,:),'b','LineWidth',2)
plot(meanHi(2,:),'b','LineWidth',1)
plot(meanHi(3,:),'b','LineWidth',1)
plot(meanLow(1,:),'r','LineWidth',2)
plot(meanLow(2,:),'r','LineWidth',1)
plot(meanLow(3,:),'r','LineWidth',1)
plot(meanRew(1,:),'b--','LineWidth',2)
plot(meanRew(2,:),'b--','LineWidth',1)
plot(meanRew(3,:),'b--','LineWidth',1)
plot(meanCatch(1,:),'k','LineWidth',2)
plot(meanCatch(2,:),'k','LineWidth',1)
plot(meanCatch(3,:),'k','LineWidth',1)
plot(meanFree(1,:),'g','LineWidth',2)
plot(meanFree(2,:),'g','LineWidth',1)
plot(meanFree(3,:),'g','LineWidth',1)
% plot([photoZero photoZero],[ylim(1) ylim(2)],'k')

set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
xlim([rasterAxis(1,2),rasterAxis(end,2)])
title('Average of Hi (b) vs Low (r) vs Aligned to Rew (k)')


%plot velocity aligned to tone
subplot(4,3,7)
hold on
plot(velRasterAxis,s.VelRaster.MeanHi(1,:),'b','LineWidth',2)
plot(velRasterAxis,s.VelRaster.MeanHi(2,:),'b','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanHi(3,:),'b','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanLow(1,:),'r','LineWidth',2)
plot(velRasterAxis,s.VelRaster.MeanLow(2,:),'r','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanLow(3,:),'r','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanFree(1,:),'g','LineWidth',2)
plot(velRasterAxis,s.VelRaster.MeanFree(2,:),'g','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanFree(3,:),'g','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanCatch(1,:),'k','LineWidth',2)
plot(velRasterAxis,s.VelRaster.MeanCatch(2,:),'k','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanCatch(3,:),'k','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanRew(1,:),'b--','LineWidth',2)
plot(velRasterAxis,s.VelRaster.MeanRew(2,:),'b--','LineWidth',1)
plot(velRasterAxis,s.VelRaster.MeanRew(3,:),'b--','LineWidth',1)
% plot([0 0],[min(min(velHiMean,velLowMean)) max(max(velHiMean,velLowMean))],'k')
title('Velocity Relative to Tone/Reward')

%Plot licking aligned to tone
subplot(4,3,10)
hold on
plot(histLickAxis,s.Licking.ToneHistHi,'b','LineWidth',2)
plot(histLickAxis,s.Licking.ToneHistLow,'r','LineWidth',2)
plot(histLickAxis,s.Licking.ToneHistCatch,'k','LineWidth',2)
plot(histLickAxis,s.Licking.ToneHistFree,'g','LineWidth',2)
title('Licking Relative to Tone')

%Plot heatmaps of photometry response
subplot(2,3,2)
imagesc(photoRaster')
colormap('parula')
set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
title(fileName)


%Plot licking rasters
subplot(4,3,8)
plot(lickRasterRew(:,1),lickRasterRew(:,2),'k.')
title('Lick Raster To Reward (0)')
xlim(rasterWindow)

subplot(4,3,11)
plot(lickRasterToneHi(:,1),lickRasterToneHi(:,2),'b.');
hold on
plot(lickRasterToneLow(:,1),lickRasterToneLow(:,2),'r.');
plot(lickRasterToneCatch(:,1),lickRasterToneCatch(:,2),'k.');
plot(lickRasterToneFree(:,1),lickRasterToneFree(:,2),'g.');

xlim(rasterWindow)
title('Lick Rasters to Tone')

% %plot out split up photometry response
% subplot(4,3,3)
% imagesc(photoRaster(:,trialLow)')
% colormap('parula')
% set(gca,'XTick',rasterAxis(:,2));
% set(gca,'XTickLabel',rasterAxis(:,1));
% title('Low Trials')
% 
% subplot(4,3,6)
% imagesc(photoRaster(:,trialHi)')
% colormap('parula')
% set(gca,'XTick',rasterAxis(:,2));
% set(gca,'XTickLabel',rasterAxis(:,1));
% title('High Trials'

%plot out photoemtry over time as binned
subplot(4,3,3)
hold on
plot(s.Vals(3,:),'b.')
plot(trialLow,s.Vals(3,trialLow),'r.')
plot(trialFree,s.Vals(3,trialFree),'g.')
plot(trialCatch,s.Vals(3,trialCatch),'k.')
xlim([0 length(toneOnset)])
title('Peak Values Across Session: hi(b) low(r)')

%plot out licking over time. 
subplot(4,3,6)
hold on
plot(binLick,'b.')
plot(trialLow,binLick(trialLow),'r.')
plot(trialFree,binLick(trialFree),'g.')
plot(trialCatch,binLick(trialCatch),'k.')
xlim([0 length(toneOnset)])
title('AntiLicks Values Across Session: hi(b) low(r)')


%plot out average velocity aligned to reward.
subplot(4,3,12)
plot(velRasterAxisRew,mean(velRasterRew'))
hold on
plot([0 0],[min(mean(velRasterRew')) max(mean(velRasterRew'))],'k')
xlim([rasterWindow(1) rasterWindow(2)])
title('Velocity Relative to Reward')


spikeGraphName = strcat(fileName,'Figure');

savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


saveName = strcat(fileName,'Analysis','.mat');
fname = saveName;
pname = pwd;

save(fullfile(pname,fname),'s');


end
















