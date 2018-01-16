%This is code that will read and link MBED inputs with photometry, while
%recording from my big rig. 
function [s] = analysisTwoTonePavlovPhotometry(fileName);
%% Parameters
rasterWindow = [-3,5]; %raster window in seconds
thresh = 0.01; %threshold for peak detection
locoTimeStep = 0.1;
minBar = 8000; %for clearing out false pulses from speaker card.
photoToggle = 0;

%use diary function to save logfile of analysis. this is good for
%troubleshooting. 
diaryName = strcat(fileName,'ANALYSISLOGFILE');
diary(diaryName)

disp(fileName)

%store things in a big structure
s = struct;
s.Parameters.RasterWindow = rasterWindow;
s.Parameters.PeakThreshold = thresh;
s.Parameters.LocomotionTimeStep = locoTimeStep;
%% First, lets pull the MBED stuff
[trialStates, portStates, trialParams,inputPhot,outputPhot,photoToggle] = functionMBEDtmp(fileName);

rewOut = [portStates.tStamps',portStates.outStates(:,2)];
toneIn = [portStates.tStamps',portStates.inStates(:,1)];

%eliminate duplicate values!
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

%use tone times for the time of tone delivery. 

[onsetPhot] = functionTTLrepairTTL(trialStates.cueTimeRew1,onsetPhot,6,1.2,100);
%assign trials to rewarded or unrewarded. first confirm that number of
%rewards is half of onsets. If not, throw error.
if length(rewTimes) ~= length(onsetPhot)/2
    disp(length(rewTimes))
    disp(length(onsetPhot))
    error('Incorrect Number of Rewards Relative to Trials')
else
    disp('Correct Number of ')
end

rewIndex = zeros(length(onsetPhot),1);
for i = 1:length(rewTimes)
    diffVal = onsetPhot - rewTimes(i);
    diffFinder = find(diffVal < 0,1,'last');
    rewIndex(diffFinder) = 1;
end

trialLow = find(rewIndex == 0);
trialHi = find(rewIndex == 1);

s.MBED.RewTimes = rewTimes;
s.MBED.RewDur = rewDur;
% s.MBED.Jitter = inputPhotOnset;
s.MBED.ToneDelivery = onsetPhot;
s.MBED.HiTrials = trialHi;
s.MBED.LowTrials = trialLow;
s.MBED.Raw = portStates;
%Now pull locomotor data

[locoData] = functionLocoTmp(fileName,portStates,locoTimeStep);

s.Locomotion = locoData;

%% Now lets pull the photometry inputs
[filtSig1,filtSig2,traceDF,traceTiming,t_ds,newSmoothDS,targetPeaks,data] = functionTDTtmp(fileName,0);

s.Photo.dFTrace = traceDF;
s.Photo.dFTime = traceTiming;
s.Photo.x70 = filtSig1;
s.Photo.x05 = filtSig2;
s.Photo.Raw = data.streams.Fi1r.data;
s.Photo.RawRate = data.streams.Fi1r.fs;
s.Photo.Peaks = targetPeaks;
s.Photo.Photo.x70dF = newSmoothDS;
s.Photo.Photo.x70dFTime = t_ds;

%pull jittered signal
traceJitt = data.epocs.PtE1.onset;

[traceJitt,inputPhotOnset] = functionTDT2MatJitterAlign(traceJitt,inputPhotOnset);
traceJittDiff = diff(traceJitt);
inputPhotDiff = diff(inputPhotOnset);

s.Photo.Jitter = traceJitt;
s.MBED.Jitter = inputPhotOnset;

%now lets deal with inputs. This is now designed to try and handle cases in
%which sound pulses are not delivered. 
% pull output timings
try
    traceMBED = data.epocs.PtC0.onset;
    interpTrig = 0;
catch
    disp('No TDT Tone Pulses Detected')
    traceMBED = interp1(inputPhotOnset,traceJitt,onsetPhot);
    interpTrig = 1;
end

%check alignment to MBED TTLs
[traceMBED] = functionTTLrepairTTL(onsetPhot/1000,traceMBED,0.3,1.2,100);
traceMBEDDiff = diff(traceMBED);

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
meanHi = mean(photoRaster(:,trialHi)');
meanLow = mean(photoRaster(:,trialLow)');
steHi = std(photoRaster(:,trialHi)')/sqrt(length(trialHi));
steLow = std(photoRaster(:,trialLow)')/sqrt(length(trialHi));
steRew = std(photoRasterRew')/sqrt(length(traceMBED));

photoZero = (-1*rasterWindow(1))/photoTimeStep;

s.PhotoRaster.ToneRaster = photoRaster;
s.PhotoRaster.RewardRaster = photoRasterRew;
s.PhotoRaster.MeanHi = meanHi;
s.PhotoRaster.MeanLow = meanLow;
s.PhotoRaster.LickRaster = photoRasterLick;
s.PhotoRaster.TimeStep = photoTimeStep;


%% Velocity Data
% Here, we will use the photometry as the base, since it is probably more
% reliable
if length(locoData.Velocity) > 0
    velTrueTime = interp1(inputPhotOnset/1000,traceJitt,locoData.Velocity(:,1));

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

s.VelRaster.ToneRaster = velRaster;
s.VelRaster.RewardRaster = velRasterRew;
s.VelRaster.Axis = velRasterAxis;

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
    lickRasterToneHi = [0,0];
    
    for i = 1:length(trialHi)
        %see if there are values that are present
        lickFinder = find(lickRasterTone(:,2) == trialHi(i));
        if length(lickFinder) > 0
            lickRasterToneHi(lickInd:lickInd + length(lickFinder)-1,:) = lickRasterTone(lickFinder,:);
            lickInd = lickInd + length(lickFinder);
        end
    end
    
    lickRasterToneLow = [0,0];

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
    lickRasterTone = zeros(1,2);
    lickRasterRew = zeros(1,2);
    lickRasterToneHi = zeros(1,2);
    lickRasterToneLow = zeros(1,2);
    histLickToneHi = zeros(length([rasterWindow(1):0.1:rasterWindow(2)]),1);
    histLickToneLow = zeros(length([rasterWindow(1):0.1:rasterWindow(2)]),1);
    histLickRewHi = zeros(length([rasterWindow(1):0.1:rasterWindow(2)]),1);
    histLickAxis = [rasterWindow(1):0.1:rasterWindow(2)];
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


s.Licking.ToneRaster = lickRasterTone;
s.Licking.RewardRaster = lickRasterRew;
s.Licking.ToneHistHi = histLickToneHi;
s.Licking.ToneHistLow = histLickToneLow;
s.Licking.ToneHistRew = histLickRewHi;
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
plot(meanHi,'b','LineWidth',2)
plot(meanHi+steHi,'b','LineWidth',1)
plot(meanHi-steHi,'b','LineWidth',1)
plot(meanLow,'r','LineWidth',2)
plot(meanLow+steLow,'r','LineWidth',1)
plot(meanLow-steLow,'r','LineWidth',1)
meanRew = mean(photoRasterRew');
plot(meanRew,'k','LineWidth',2)
plot(meanRew + steRew,'k')
plot(meanRew - steRew,'k')
plot([photoZero photoZero],[min(min(meanHi,meanLow)) max(max(meanHi,meanLow))],'k')

set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
xlim([rasterAxis(1,2),rasterAxis(end,2)])
title('Average of Hi (b) vs Low (r) vs Aligned to Rew (k)')




%plot velocity aligned to tone
subplot(4,3,7)
hold on
plot(velRasterAxis,velHiMean,'b','LineWidth',2)
plot(velRasterAxis,velHiMean+velHiSTE,'b','LineWidth',1)
plot(velRasterAxis,velHiMean-velHiSTE,'b','LineWidth',1)
plot(velRasterAxis,velLowMean,'r','LineWidth',2)
plot(velRasterAxis,velLowMean+velLowSTE,'r','LineWidth',1)
plot(velRasterAxis,velLowMean-velLowSTE,'r','LineWidth',1)

plot([0 0],[min(min(velHiMean,velLowMean)) max(max(velHiMean,velLowMean))],'k')
title('Velocity Relative to Tone')

%Plot licking aligned to tone
subplot(4,3,10)
hold on
plot(histLickAxis,histLickToneHi,'b','LineWidth',2)
plot(histLickAxis,histLickToneLow,'r','LineWidth',2)
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
plot(lickRasterToneLow(:,1),lickRasterToneLow(:,2),'r.')

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
title('Peak Values Across Session: hi(b) low(r)')

%plot out licking over time. 
subplot(4,3,6)
hold on
plot(binLick,'b.')
plot(trialLow,binLick(trialLow),'r.')
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
diary off

end

















