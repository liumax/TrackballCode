%This is code that will read and link MBED inputs with photometry, while
%recording from my big rig. 

%% Parameters

rasterWindow = [-3,5]; %raster window in seconds
thresh = 0.01;
%% First, lets pull the MBED stuff

% 
% [portStates] = maxTrialVariableNoTask(fileName);
% 
% %we want to look at instates. Port 1 is TDT for photometry, port 2 is
% %NOLDUS
% 
% inputPhot = [portStates.tStamps',portStates.inStates(:,8)];
% outputPhot = [portStates.tStamps',portStates.outStates(:,8)];
% 
% %now eliminate duplicates in photometry inputs
% whileTrig = 1;
% whileInd = 2;
% 
% while whileTrig == 1;
%     if whileInd > length(inputPhot)
%         break
%     end
%     prevVal = inputPhot(whileInd-1,2);
%     currVal = inputPhot(whileInd,2);
%     if prevVal == currVal
% %         disp('Duplicate Detected!')
%         inputPhot(whileInd,:) = [];
%     else
%         whileInd = whileInd + 1;
%     end
% end
% 
% inputPhotOnset = inputPhot(inputPhot(:,2) == 1,1);
% 
% %clean up input signal to remove excess. 
% inputDiff = diff(inputPhotOnset);
% diffFind = find(inputDiff > 1000);
% if length(diffFind) == 1
%     inputPhotOnset(1:diffFind) = [];
% else
%     error('multiple input time periods')
% end
% 
% 
% %do the same for photometry outputs
% whileTrig = 1;
% whileInd = 2;
% 
% while whileTrig == 1;
%     if whileInd > length(outputPhot)
%         break
%     end
%     prevVal = outputPhot(whileInd-1,2);
%     currVal = outputPhot(whileInd,2);
%     if prevVal == currVal
% %         disp('Duplicate Detected!')
%         outputPhot(whileInd,:) = [];
%     else
%         whileInd = whileInd + 1;
%     end
% end
% 
% %pull just the onset times for photometry
% onsetPhot = outputPhot(outputPhot(:,2) == 1,1);
% onsetPhotDiff = diff(onsetPhot);
% 
% %Now pull locomotor data
% 
% [locoData] = functionMBEDrotary(portStates.inStates(:,4),portStates.inStates(:,5),portStates.tStamps/1000,0.1);

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

% [peakInfo, riseInfo, troughInfo] = findPhotoPeaks(traceTiming,traceDF,thresh);

%pull jittered signal
traceJitt = data.epocs.PtE1.onset;
traceJittDiff = diff(traceJitt);
% 
% %fix the MBED jittered signal
% if length(inputPhotOnset) > length(traceJitt);
%     disp('More MBED INPUTS IN JITTER, SUBTRACTING')
%     inputPhotOnset(length(traceJitt)+1:end) = [];
% elseif length(inputPhotOnset) < length(traceJitt);
%     error('Fewer MBED INPUTS THAN REPORTED')
% elseif length(inputPhotOnset) == length(traceJitt);
%     disp('Jittered traces matched!')
% end
% 
%pull output timings
try
    traceMBED = data.epocs.PtC0.onset;
catch
    disp('No TDT Tone Pulses Detected')
%     traceMBED = interp1(inputPhotOnset,traceJitt,onsetPhot);
end
traceMBEDDiff = diff(traceMBED);
% 
% %check alignment
% [xcf,lags,bounds]  = crosscorr(onsetPhotDiff,traceMBEDDiff);
% [xcMax maxInd] = max(xcf);
% xcLag = lags(maxInd);
% 
% if xcLag ~= 0
%     error('Photometry Not Aligned')
% elseif length(onsetPhot) ~= length(traceMBED)
%     error('Mismatch in Number of Photometry Pulses')
% end

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
% [riseRasters] = functionBasicRaster((riseInfo.t)',traceMBED,rasterWindow);
%generate correct order for displaying things by freq/db
sortingCounter = 1;
for i = 1:numFreqs
    for j = 1:numDBs
        sortingFinder = find(trialMatrix(:,2) == uniqueFreqs(i) & trialMatrix(:,3) == uniqueDBs(j));
        sortIndex(sortingFinder) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end

% riseRasters(:,3) = 0;
% riseRasters(:,3) = sortIndex(riseRasters(:,2));

%% Interpolate times for cross plotting of MBED data and Photometry Signal
% Here, we will use the photometry as the base, since it is probably more
% reliable

% velTrueTime = interp1(onsetPhot/1000,traceMBED,locoData.Velocity(:,1));
% 
% %in the current iteration, this has issues because the photometry cant
% %align without pulses, and i only have output pulses from the MBED when the
% %sound turns on. Therefore, I will need to crop things out.
% 
% findVelFirst = find(~isnan(velTrueTime),1,'first');
% findVelLast = find(~isnan(velTrueTime),1,'last');
% %find nearest in photometry signal
% findPhotFirst = find(traceTiming - velTrueTime(findVelFirst)>0,1,'first');
% findPhotLast = find(traceTiming - velTrueTime(findVelLast)>0,1,'first');

        
%% Plot everything
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
%create a set of ticks for labeling any display of octaves
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

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])

%plot overall photometry trace and locomotion trace
% subplot(3,3,1)
% hold on
% % plot(velTrueTime(findVelFirst:findVelLast),(locoData.Velocity(findVelFirst:findVelLast,2)-min(locoData.Velocity(findVelFirst:findVelLast,2)))/(max(locoData.Velocity(findVelFirst:findVelLast,2))-min(locoData.Velocity(findVelFirst:findVelLast,2))))
% plot(traceTiming(findPhotFirst:1000:findPhotLast),(traceDF(findPhotFirst:1000:findPhotLast)-min(traceDF(findPhotFirst:findPhotLast)))/(max(traceDF(findPhotFirst:findPhotLast))-min(traceDF(findPhotFirst:findPhotLast))),'r')
% xlim([velTrueTime(findVelFirst),velTrueTime(findVelLast)])
% title('Normalized Vel (b) and Photometry (r)')

%plot rasters?

%Plot heatmaps of photometry response
subplot(3,3,3)
imagesc(squeeze(photoAverages(:,:,1)'))
colormap('parula')
set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
set(gca,'YTick',octaveRange(:,2));
set(gca,'YTickLabel',octaveRange(:,1));
title({fileName;'60DB'})
        
subplot(3,3,6)
imagesc(squeeze(photoAverages(:,:,2)'))
colormap('parula')
set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
set(gca,'YTick',octaveRange(:,2));
set(gca,'YTickLabel',octaveRange(:,1));
title('80DB')

subplot(3,3,9)
imagesc(squeeze(photoAverages(:,:,3)'))
colormap('parula')
set(gca,'XTick',rasterAxis(:,2));
set(gca,'XTickLabel',rasterAxis(:,1));
set(gca,'YTick',octaveRange(:,2));
set(gca,'YTickLabel',octaveRange(:,1));
title('100DB')

spikeGraphName = strcat(fileName,'Figure');

savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
























