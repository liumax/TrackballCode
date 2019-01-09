

rasterWindow = [-3 5];

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
        testDiff = abs(rewTimes - onsetPhot(i))
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
end


lickData = trialParams.licking;


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


figure
%Plot licking aligned to tone
subplot(3,1,1)
hold on
plot(histLickAxis,histLickToneHi,'b','LineWidth',2)
plot(histLickAxis,histLickToneLow,'r','LineWidth',2)
title('Licking Relative to Tone')

%Plot licking rasters
subplot(3,1,2)
plot(lickRasterRew(:,1),lickRasterRew(:,2),'k.')
title('Lick Raster To Reward (0)')
xlim(rasterWindow)

subplot(3,1,3)
plot(lickRasterToneHi(:,1),lickRasterToneHi(:,2),'b.');
hold on
plot(lickRasterToneLow(:,1),lickRasterToneLow(:,2),'r.')