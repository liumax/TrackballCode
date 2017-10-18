%This is code that will read and link MBED inputs with photometry, while
%recording from my big rig. 

function [s] = analysisPhotometryTuning(fileName);

%% Parameters

rasterWindow = [-3,5]; %raster window in seconds
thresh = 0.01;
locoTimeStep = 0.1;
photoToggle = 0;


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
try
    [inputPhot] = functionSignalDuplicateElim(inputPhot,2);
catch
    disp('NO INPUTS, SWITCH TO USING PHOTOMETRY')
    photoToggle = 1;
end
[outputPhot] = functionSignalDuplicateElim(outputPhot,2);


%pull just the onset times for photometry
onsetPhot = outputPhot(outputPhot(:,2) == 1,1);
onsetPhotDiff = diff(onsetPhot);


try
    inputPhotOnset = inputPhot(inputPhot(:,2) == 1,1);
    inputPhotDiff = diff(inputPhotOnset);
catch
    disp('NO INPUTS, CANT PULL ONSETS')
end


s.MBED.ToneDelivery = onsetPhot;

%Now pull locomotor data

[locoData] = functionMBEDrotary(portStates.inStates(:,4),portStates.inStates(:,5),portStates.tStamps/1000,locoTimeStep);

s.Locomotion = locoData;

%% Now lets pull the photometry inputs
%load file
data = load(strcat(fileName,'.mat'));
data=data.data;

[filtSig1,filtSig2,traceDF,traceTiming] = functionPhotometryRawExtraction(data);

s.Photo.dFTrace = traceDF;
s.Photo.dFTime = traceTiming;
s.Photo.x70 = filtSig1;
s.Photo.x05 = filtSig2;
% s.Photo.Raw = data.streams.Fi1r.data;
% s.Photo.RawRate = data.streams.Fi1r.fs;

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

if photoToggle == 1
    disp('Generating Fake Jitter Based on MBED Tone Pulses')
    traceMBED = data.epocs.PtC0.onset;
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
    
    %generate input times off of interpolation
    inputPhotOnset = interp1(traceMBED,onsetPhot,traceJitt);
    disp('Fake Jitter Constructed...')
    tester = find(isnan(inputPhotOnset));
    traceJitt(tester) = [];
    inputPhotOnset(tester) = [];
    disp('Deleting NaN values')
    traceJittDiff = diff(traceJitt);
    inputPhotDiff = diff(inputPhotOnset);
end


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
if length(massDiff) <=10 & length(massDiff)>0
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
end

s.Photo.Jitter = traceJitt;
s.MBED.Jitter = inputPhotOnset;


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

s.Photo.MBEDSig = traceMBED; %store tone times, makes life easier later on. 

%calculate raster in terms of time steps in photometry
photoTimeStep = mean(diff(t_ds));
rasterPhotWindow = round(rasterWindow/photoTimeStep);

s.Photo.AlignTimes = traceMBED;
                                                        
%% Pull sound data!

%extract data into matlab
soundName = strcat(fileName,'Sound.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info

s.SoundData = soundData;

%I basically want to just figure out the order of things, since I should
%have the pulses inside the photometry trace already. This should be found
%in trial matrix. 

trialMatrix = soundData.TrialMatrix;

uniqueFreqs = unique(trialMatrix(:,2));
uniqueDBs = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(1),3));

numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);

%now lets also sort out whether this is the right tuning curve. do
%crosscorr on the delays vs outputs. 

[xcf,lags,bounds]  = crosscorr(onsetPhotDiff,soundData.Delays);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);

if xcLag == 0 & xcMax > 0.90
    disp('Sound File and Actual Lags Aligned!')
elseif xcLag < 0 & xcMax > 0.9
    disp('Negative Lag detected in Sound File vs Real Data. Deleting data from OnsetPhot')
    onsetPhot(1:-xcLag) = [];
    onsetPhotDiff = diff(onsetPhot);
    traceMBED(1:-xcLag) = [];
    s.Photo.MBEDSig = traceMBED; %store tone times, makes life easier later on. 
    s.Photo.AlignTimes = traceMBED;

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

for ind = 1:length(traceMBED)
    alignTime = traceMBED(ind);
    %find the time in the photometry trace
    photoPoint = find(t_ds - alignTime > 0,1,'first');
    if photoPoint + rasterPhotWindow(2) < length(newSmoothDS)
        photoRaster(:,ind) = newSmoothDS(photoPoint + rasterPhotWindow(1):photoPoint + rasterPhotWindow(2));
    else
        disp('Tone Rasters: Reached End of Photometry Trace')
        disp(ind)
        break
    end
    
end

%make simple plot
simplePlot = photoRaster;
for ind = 1:length(traceMBED)
    simplePlot(:,ind) = simplePlot(:,ind)+ind*0.5;
end

s.Processed.PhotoRaster = photoRaster;


%now I need to store different means

photoAverages = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,numFreqs,numDBs);
photoRasterStore = cell(numFreqs,numDBs);
for ind = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(ind),3));
    for j = 1:numDBs
        %find all trials of the particular setting
        targetFinder = find(trialMatrix(:,2) == uniqueFreqs(ind) & trialMatrix(:,3) == subUniqueDB(j));
        %pull and average these traces
        tempHolder = photoRaster(:,targetFinder);
        photoRasterStore{ind,j} = tempHolder;
        tempHolder = mean(tempHolder');
        photoAverages(:,ind,j) = tempHolder;
    end
end

s.Processed.PhotoAverages = photoAverages;


%lets make rasters with the peak detection rise times (based on Chris'
%recommendation)

%make basic rasters
[riseRasters] = functionBasicRaster((targetPeaks(:,4)),traceMBED,rasterWindow);
%generate correct order for displaying things by freq/db
sortingCounter = 1;
for ind = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(ind),3));
    for j = 1:numDBs
        sortingFinder = find(trialMatrix(:,2) == uniqueFreqs(ind) & trialMatrix(:,3) == subUniqueDB(j));
        sortIndex(sortingFinder) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end
try
    riseRasters(:,3) = 0;
    riseRasters(:,3) = sortIndex(riseRasters(:,2));
catch
    riseRasters(:,3) = 0;
end

s.Processed.RiseRaster = riseRasters;

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

%find nearest in photometry signal
velTrueTime(isnan(velTrueTime)) = [];
if velTrueTime
    findVelFirst = find(~isnan(velTrueTime),1,'first');
    findVelLast = find(~isnan(velTrueTime),1,'last');
    findPhotFirst = find(traceTiming - velTrueTime(findVelFirst)>0,1,'first');
    findPhotLast = find(traceTiming - velTrueTime(findVelLast)>0,1,'first');
else
    findVelFirst = 1;
    findVelLast = length(traceTiming);
    findPhotFirst = 1;
    findPhotLast = length(traceTiming);
end

        
%% Plot everything
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

% create a set of ticks for labeling any display of octaves
if soundData.WhiteNoise
    if soundData.WhiteNoise == 1
        totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(2)));
        
        %This then makes an array of the full octave steps I've made
        octaveRange = zeros(totalOctaves + 1,2);
        octaveRange(1,1) = uniqueFreqs(2);
        for i = 1:totalOctaves
            octaveRange (i+1,1) = octaveRange(i,1)*2;
        end
        %next, I find the positions from uniqueFreqs that match octaveRange
        for i = 1:size(octaveRange,1);
            octaveRange(i,2) = find(round(uniqueFreqs) == octaveRange(i,1));
        end
        %add in the white noise
        octaveRange(2:end+1,:) = octaveRange(1:end,:);
        octaveRange(1,1) = 0;
        octaveRange(1,2) = 1;
    else soundData.WhiteNoise == 0;
        totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));
        
        %This then makes an array of the full octave steps I've made
        octaveRange = zeros(totalOctaves + 1,2);
        octaveRange(1,1) = uniqueFreqs(1);
        for ind = 1:totalOctaves
            octaveRange (ind+1,1) = octaveRange(ind,1)*2;
        end
        %next, I find the positions from uniqueFreqs that match octaveRange
        for ind = 1:size(octaveRange,1);
            octaveRange(ind,2) = find(round(uniqueFreqs) == octaveRange(ind,1));
        end
    end
else
    totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));
    
    %This then makes an array of the full octave steps I've made
    octaveRange = zeros(totalOctaves + 1,2);
    octaveRange(1,1) = uniqueFreqs(1);
    for ind = 1:totalOctaves
        octaveRange (ind+1,1) = octaveRange(ind,1)*2;
    end
    %next, I find the positions from uniqueFreqs that match octaveRange
    for ind = 1:size(octaveRange,1);
        octaveRange(ind,2) = find(round(uniqueFreqs) == octaveRange(ind,1));
    end
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
for ind = 1:numFreqs
    rasterTicks{ind} = num2str(rasterLabels(ind,1));
end

%MAKE THE FIGURE
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])

%plot overall photometry trace and locomotion trace
subplot(3,3,1)
hold on
if velTrueTime
    plot(velTrueTime(findVelFirst:findVelLast),(locoData.Velocity(findVelFirst:findVelLast,2)-min(locoData.Velocity(findVelFirst:findVelLast,2)))/(max(locoData.Velocity(findVelFirst:findVelLast,2))-min(locoData.Velocity(findVelFirst:findVelLast,2))))
    xlim([velTrueTime(findVelFirst),velTrueTime(findVelLast)])
end
plot(traceTiming(findPhotFirst:1000:findPhotLast),(traceDF(findPhotFirst:1000:findPhotLast)-min(traceDF(findPhotFirst:findPhotLast)))/(max(traceDF(findPhotFirst:findPhotLast))-min(traceDF(findPhotFirst:findPhotLast))),'r')

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
for ind = 1:numFreqs
    plot([rasterWindow(1) rasterWindow(2)],[rasterLabels(ind,2) rasterLabels(ind,2)],'g')
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

% %now make a figure of individual traces!
% 
% %first, determine max and min values overall
% 
% for i = 1:numFreqs
%     for j = 1:numDBs
%         openFile = photoRasterStore{i,j};
%         minVal(i,j) = min(min(openFile));
%         maxVal(i,j) = max(max(openFile));
%     end
% end
% 
% fullMinVal = min(min(minVal));
% fullMaxVal = max(max(maxVal));
% 
% %find the zero point on the raster
% zeroPoint = find(rasterAxis(:,1) == 0);
% zeroPoint = rasterAxis(zeroPoint,2);
% 
% hFig = figure;
% 
% for i = 1:numDBs
%     subplot(1,numDBs,i)
%     hold on
%     for j = 1:numFreqs
%         plot(photoRasterStore{j,i} + fullMaxVal*j)
%         bigMax = max(max(photoRasterStore{j,i} + fullMaxVal*j));
%     end
%     xlim([0 length(photoRasterStore{j,i})])
%     set(gca,'XTick',rasterAxis(:,2));
%     set(gca,'XTickLabel',rasterAxis(:,1));
%     plot([zeroPoint zeroPoint],[fullMinVal+fullMaxVal bigMax],'LineWidth',2)
% end


saveName = strcat(fileName,'Analysis','.mat');
fname = saveName;
pname = pwd;

save(fullfile(pname,fname),'s');














end


