%This is code that will read and link MBED inputs with photometry, while
%recording from my big rig. 

function [s] = analysisPhotometryTuning(fileName);

%% Parameters

rasterWindow = [-3,5]; %raster window in seconds
thresh = 0.01;
locoTimeStep = 0.1;


%store things in a big structure
s = struct;
s.Parameters.RasterWindow = rasterWindow;
s.Parameters.PeakThreshold = thresh;
s.Parameters.LocomotionTimeStep = locoTimeStep;


%% First, lets pull the MBED stuff


[trialStates, portStates, trialParams] = maxTrialVariablesLickingTask(fileName);

%we want to look at instates. Port 1 is TDT for photometry, port 2 is
%NOLDUS

inputPhot = [portStates.tStamps',portStates.inStates(:,8)];
outputPhot = [portStates.tStamps',portStates.outStates(:,8)];

%eliminate duplicate values!
[inputPhot] = functionSignalDuplicateElim(inputPhot,2);
[outputPhot] = functionSignalDuplicateElim(outputPhot,2);


%pull just the onset times for photometry
onsetPhot = outputPhot(outputPhot(:,2) == 1,1);
onsetPhotDiff = diff(onsetPhot);

inputPhotOnset = inputPhot(inputPhot(:,2) == 1,1);
inputPhotDiff = diff(inputPhotOnset);

s.MBED.ToneDelivery = onsetPhot;

%Now pull locomotor data

[locoData] = functionMBEDrotary(portStates.inStates(:,4),portStates.inStates(:,5),portStates.tStamps/1000,locoTimeStep);

s.Locomotion = locoData;

%% Now lets pull the photometry inputs
%load file
data = load(strcat(fileName,'.mat'));
data=data.data;


%Code from Chris that performs isosbestic correction. 
traceDF = isosbestic_correction(data);


%pull timestamps for fluorescence
traceTiming = [0:1/data.streams.x70G.fs:(1/data.streams.x70G.fs)*(length(data.streams.x70G.data)-1)];

s.Photo.dFTrace = traceDF;
s.Photo.dFTime = traceTiming;
%pull peaks 170616 This appears to have problem: built for 2016 matlab, has
%additional functionality for peak finding.
try
    [peakInfo, riseInfo, troughInfo] = findPhotoPeaks(traceTiming,traceDF,thresh);
catch
    disp('FIND PEAKS CODE FAILED')
    peakInfo = [];
    riseInfo = [];
    troughInfo = [];
end

s.Photo.Peaks.Peak = peakInfo;
s.Photo.Peaks.Rise = riseInfo;
s.Photo.Peaks.Trough = troughInfo;


%pull jittered signal
traceJitt = data.epocs.PtE1.onset;
traceJittDiff = diff(traceJitt);


%clean up the MBED signal. generally there will be excess TTLs at the
%beginning
findBig = find(inputPhotDiff > 600);
%now we need to screen the big differences.
for i = 1:length(findBig)
    bigSize = findBig(i);
    if bigSize > length(inputPhotOnset)/2 %in the case of something coming near the end
        disp('Late Big Diff, deleting')
        inputPhotOnset(bigSize:end) = [];
        inputPhotDiff = diff(inputPhotOnset);
        findBig = find(inputPhotDiff > 600);
    else
        disp('Early Big Difference in Jitter')
    end
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
[xcf,lags,bounds]  = crosscorr(inputPhotDiff,traceJittDiff);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);

if xcLag ~= 0
    disp('CorrLag')
    disp(xcLag)
    disp('MaxCorr')
    disp(xcMax)
    error('Jitter Not Aligned')
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
elseif length(inputPhotDiff) == length(traceJittDiff)
    disp('TRACES ARE ALIGNED!')
end

s.Photo.Jitter = traceJitt;
s.MBED.Jitter = inputPhot;


%now lets deal with inputs. This is now designed to try and handle cases in
%which sound pulses are not delivered. 
% pull output timings
try
    traceMBED = data.epocs.PtC0.onset;
catch
    disp('No TDT Tone Pulses Detected')
    traceMBED = interp1(inputPhotOnset,traceJitt,onsetPhot);
end
traceMBEDDiff = diff(traceMBED);

%check alignment
[xcf,lags,bounds]  = crosscorr(onsetPhotDiff/1000,traceMBEDDiff);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);

if xcLag ~= 0
    error('Tones Not Aligned')
elseif length(onsetPhot) ~= length(traceMBED)
    error('Mismatch in Number of Tone Pulses')
end

%calculate raster in terms of time steps in photometry
photoTimeStep = 1/data.streams.x70G.fs;
rasterPhotWindow = round(rasterWindow/photoTimeStep);

%% Pull sound data!

%extract data into matlab
soundName = strcat(fileName,'Sound.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info


%I basically want to just figure out the order of things, since I should
%have the pulses inside the photometry trace already. This should be found
%in trial matrix. 

trialMatrix = soundData.TrialMatrix;

uniqueFreqs = unique(trialMatrix(:,2));
uniqueDBs = unique(trialMatrix(:,3));

numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);

%now lets also sort out whether this is the right tuning curve. do
%crosscorr on the delays vs outputs. 

[xcf,lags,bounds]  = crosscorr(onsetPhotDiff,soundData.Delays);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);

if xcLag == 0 & xcMax > 0.90
    disp('Sound File and Actual Lags Aligned!')
else
    error('Check if Correct Tuning File was Used!')
end

%% Now lets sort out tuning curve

%First, check to see if number of tones equals number of outputs.

if length(trialMatrix) == length(traceMBED)
    disp('Trial Matrix and Photometry Inputs Matched')
elseif length(trialMatrix) ~= length(traceMBED)
    error('Mismatch in Trial Number and Photometry Inputs')
end

numReps = soundData.ToneRepetitions;
toneDurPhot = round(soundData.ToneDuration/photoTimeStep);

%first, lets just pull the basic raster

photoRaster = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,length(traceMBED));

for i = 1:length(traceMBED)
    alignTime = traceMBED(i);
    %find the time in the photometry trace
    photoPoint = find(traceTiming - alignTime > 0,1,'first');
    photoRaster(:,i) = traceDF(photoPoint + rasterPhotWindow(1):photoPoint + rasterPhotWindow(2));
end

%make simple plot
simplePlot = photoRaster;
for i = 1:length(traceMBED)
    simplePlot(:,i) = simplePlot(:,i)+i*0.5;
end

%now I need to store different means

photoAverages = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,numFreqs,numDBs);
photoRasterStore = cell(numFreqs,numDBs);
for i = 1:numFreqs
    for j = 1:numDBs
        %find all trials of the particular setting
        targetFinder = find(trialMatrix(:,2) == uniqueFreqs(i) & trialMatrix(:,3) == uniqueDBs(j));
        %pull and average these traces
        tempHolder = photoRaster(:,targetFinder);
        photoRasterStore{i,j} = tempHolder;
        tempHolder = mean(tempHolder');
        photoAverages(:,i,j) = tempHolder;
    end
end

%lets make rasters with the peak detection rise times (based on Chris'
%recommendation)

%make basic rasters
[riseRasters] = functionBasicRaster((riseInfo.t)',traceMBED,rasterWindow);
%generate correct order for displaying things by freq/db
sortingCounter = 1;
for i = 1:numFreqs
    for j = 1:numDBs
        sortingFinder = find(trialMatrix(:,2) == uniqueFreqs(i) & trialMatrix(:,3) == uniqueDBs(j));
        sortIndex(sortingFinder) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end

riseRasters(:,3) = 0;
riseRasters(:,3) = sortIndex(riseRasters(:,2));

% %now lets try and make a gaussian convolution of riseRasters!
% gaussAverage = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,numFreqs,numDBs);
% gaussStore = cell(numFreqs,numDBs);
% for i = 1:numFreqs
%     for j = 1:numDBs
%         %find all trials of the particular setting
%         targetFinder = find(trialMatrix(:,2) == uniqueFreqs(i) & trialMatrix(:,3) == uniqueDBs(j));
%         %pull and average these traces
%         tempCount = 1;
%         for k = 1:length(targetFinder)
%             targetFind = find(riseRasters(:,2) == targetFinder(k));
%             tempHold(tempCount:tempCount+length(targetFind)-1) = targetFind;
%             tempCount = tempCount+length(targetFind);
%         end
%         gaussStore{i,j} = riseRasters(tempHold,:);
%         %now collapse to a single line
%         newGauss = sort(riseRasters(tempHold,1));
%         photoAverages(:,i,j) = tempHolder;
%     end
% end


%% Interpolate times for cross plotting of MBED data and Photometry Signal
% Here, we will use the photometry as the base, since it is probably more
% reliable

velTrueTime = interp1(onsetPhot/1000,traceMBED,locoData.Velocity(:,1));
% 
% %in the current iteration, this has issues because the photometry cant
% %align without pulses, and i only have output pulses from the MBED when the
% %sound turns on. Therefore, I will need to crop things out.
% 
findVelFirst = find(~isnan(velTrueTime),1,'first');
findVelLast = find(~isnan(velTrueTime),1,'last');
%find nearest in photometry signal
findPhotFirst = find(traceTiming - velTrueTime(findVelFirst)>0,1,'first');
findPhotLast = find(traceTiming - velTrueTime(findVelLast)>0,1,'first');

        
%% Plot everything
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

% create a set of ticks for labeling any display of octaves
totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));
%This then makes an array of the full octave steps I've made
octaveRange = zeros(totalOctaves + 1,2);
octaveRange(1,1) = uniqueFreqs(1);
for i = 1:totalOctaves
    octaveRange (i+1,1) = octaveRange(i,1)*2;
end
%next, I find the positions from uniqueFreqs that match octaveRange
for i = 1:size(octaveRange,1);
    octaveRange(i,2) = find(round(uniqueFreqs) == octaveRange(i,1));
end
%create a set of points for the raster axis
rasterAxis = zeros(rasterWindow(2)-rasterWindow(1)+1,2);
rasterAxis(:,1) = [rasterWindow(1):rasterWindow(2)];
%find a second in photometry sample time
rasterAxis(:,2) = [0:length(photoRaster)/8:length(photoRaster)];
rasterAxis(1,2) = 1;

%also create labels for the rasters. first column is frequencies, second is
%where lines are drawn, third is for placement of labels
rasterLabels = zeros(numFreqs,3);
rasterLabels(:,1) = uniqueFreqs;
rasterLabels(:,2) = [numDBs*soundData.ToneRepetitions:numDBs*soundData.ToneRepetitions:length(soundData.TrialMatrix)];
rasterLabels(:,3) = [numDBs*soundData.ToneRepetitions/2:numDBs*soundData.ToneRepetitions:length(soundData.TrialMatrix)-numDBs*soundData.ToneRepetitions/2];

rasterTicks = cell(numFreqs,1);
for i = 1:numFreqs
    rasterTicks{i} = num2str(rasterLabels(i,1));
end

%MAKE THE FIGURE
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])

%plot overall photometry trace and locomotion trace
subplot(3,3,1)
hold on
plot(velTrueTime(findVelFirst:findVelLast),(locoData.Velocity(findVelFirst:findVelLast,2)-min(locoData.Velocity(findVelFirst:findVelLast,2)))/(max(locoData.Velocity(findVelFirst:findVelLast,2))-min(locoData.Velocity(findVelFirst:findVelLast,2))))
plot(traceTiming(findPhotFirst:1000:findPhotLast),(traceDF(findPhotFirst:1000:findPhotLast)-min(traceDF(findPhotFirst:findPhotLast)))/(max(traceDF(findPhotFirst:findPhotLast))-min(traceDF(findPhotFirst:findPhotLast))),'r')
xlim([velTrueTime(findVelFirst),velTrueTime(findVelLast)])
title('Normalized Vel (b) and Photometry (r)')

%Plot heatmaps of photometry response
%First, find the overall max and min, so that I can set the same limits for
%all imagesc
imagescLim = [min(min(min(photoAverages))),max(max(max(photoAverages)))];

subplot(3,3,2)
imagesc(squeeze(photoAverages(:,:,1)'),imagescLim)
colormap('parula')
colorbar
set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
set(gca,'YTick',octaveRange(:,2));
set(gca,'YTickLabel',octaveRange(:,1));
title({fileName;'60DB'})
        
subplot(3,3,5)
imagesc(squeeze(photoAverages(:,:,2)'),imagescLim)
colormap('parula')
colorbar
set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
set(gca,'YTick',octaveRange(:,2));
set(gca,'YTickLabel',octaveRange(:,1));
title('80DB')

subplot(3,3,8)
imagesc(squeeze(photoAverages(:,:,3)'),imagescLim)
colormap('parula')
colorbar
set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
set(gca,'YTick',octaveRange(:,2));
set(gca,'YTickLabel',octaveRange(:,1));
title('100DB')


subplot(3,3,3)
plot(riseRasters(:,1),riseRasters(:,2),'k.')
hold on
plot([0 0],[1 length(trialMatrix)],'b')
plot([soundData.ToneDuration soundData.ToneDuration],[1 length(trialMatrix)],'b')
title('Rasters Organized by Time')
ylim([1 length(trialMatrix)])
xlim(rasterWindow)

subplot(3,3,6)
plot(riseRasters(:,1),riseRasters(:,3),'k.')
hold on
for i = 1:numFreqs
    plot([rasterWindow(1) rasterWindow(2)],[rasterLabels(i,2) rasterLabels(i,2)],'g')
end
plot([0 0],[1 length(trialMatrix)],'b')
plot([soundData.ToneDuration soundData.ToneDuration],[1 length(trialMatrix)],'b')
title('Rasters Organized by Freq/DB')
ylim([1 length(trialMatrix)])
xlim(rasterWindow)
set(gca,'YTick',rasterLabels(:,3))
set(gca,'YTickLabel',rasterTicks)
set(gca,'Ydir','reverse')



spikeGraphName = strcat(fileName,'Figure');

savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
























