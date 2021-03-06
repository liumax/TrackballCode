%This is code that will read and link MBED inputs with photometry, while
%recording from my big rig. 

function [s] = analysisPhotometryTuning(fileName);

%% Parameters

rasterWindow = [-3,5]; %raster window in seconds
thresh = 0.01;
locoTimeStep = 0.1;
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
%% Next, lets pull the MBED stuff
[trialStates, portStates, trialParams,inputPhot,outputPhot,photoToggle] = functionMBEDtmp(fileName);

%pull just the onset times for photometry
onsetPhot = outputPhot(outputPhot(:,2) == 1,1);
onsetPhotDiff = diff(onsetPhot);


try
    inputPhotOnset = inputPhot(inputPhot(:,2) == 1,1);
    inputPhotDiff = diff(inputPhotOnset);
catch
    disp('NO INPUTS, CANT PULL ONSETS')
end

[onsetPhot] = functionITIrepairTTL(soundData.Delays*1000,onsetPhot,500,1.4,100);
% onsetPhot = onsetPhot * 1000;

s.MBED.ToneDelivery = onsetPhot;

%Now pull locomotor data, only if no TMP file
[locoData] = functionLocoTmp(fileName,portStates,locoTimeStep);
s.Locomotion = locoData;

%% Now lets pull the photometry inputs
%Check for tmp file!
[filtSig1,filtSig2,traceDF,traceTiming,t_ds,newSmoothDS,targetPeaks,data] = functionTDTtmp(fileName,0);

s.Photo.dFTrace = traceDF;
s.Photo.dFTime = traceTiming;
s.Photo.x70 = filtSig1;
s.Photo.x05 = filtSig2;
s.Photo.Peaks = targetPeaks;
s.Photo.Photo.x70dF = newSmoothDS;
s.Photo.Photo.x70dFTime = t_ds;

%pull jittered signal
traceJitt = data.epocs.PtE1.onset;
traceJittDiff = diff(traceJitt);

if photoToggle == 0
    [traceJitt,inputPhotOnset] = functionTDT2MatJitterAlign(traceJitt,inputPhotOnset);
else
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

[traceMBED] = functionTTLrepairTTL(onsetPhot,traceMBED*1000,20,1.4,100);
traceMBED = traceMBED/1000;
traceMBEDDiff = diff(traceMBED);

%calculate raster in terms of time steps in photometry

s.Photo.MBEDSig = traceMBED; %store tone times, makes life easier later on. 

%calculate raster in terms of time steps in photometry
photoTimeStep = mean(diff(t_ds));
rasterPhotWindow = round(rasterWindow/photoTimeStep);
rasterVelWindow = round(rasterWindow/locoTimeStep);
rasterVect = [rasterWindow(1):photoTimeStep:rasterWindow(2)];
zeroPoint = find(rasterVect > 0,1,'first');

s.Photo.AlignTimes = traceMBED;

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
velRaster = zeros(rasterVelWindow(2) - rasterVelWindow(1) + 1,length(traceMBED));
for ind = 1:length(traceMBED)
    alignTime = traceMBED(ind);
    %find the time in the photometry trace
    photoPoint = find(t_ds - alignTime > 0,1,'first');
    if photoPoint + rasterPhotWindow(2) < length(newSmoothDS)
        photoRaster(:,ind) = newSmoothDS(photoPoint + rasterPhotWindow(1):photoPoint + rasterPhotWindow(2));
    else
        disp('Tone Rasters: Reached End of Photometry Trace')
        disp(ind)
    end
    %find the time in the velocity trace
    alignTime = onsetPhot(ind)/1000;
    velPoint = find(locoData.Velocity(:,1) - alignTime > 0,1,'first');
    if velPoint + rasterVelWindow(2) < length(locoData.Velocity(:,1)) & velPoint + rasterVelWindow(1) >= 1
        velRaster(:,ind) = locoData.Velocity(velPoint + rasterVelWindow(1):velPoint + rasterVelWindow(2),2); 
    else
        disp(ind)
        disp('Vel Rasters: Not matched up')
    end
    
end

[photoRasterZ,baselineMean,baselineSTD] = functionZScore(photoRaster,zeroPoint,length(traceMBED));

%make simple plot
simplePlot = photoRaster;
for ind = 1:length(traceMBED)
    simplePlot(:,ind) = simplePlot(:,ind)+ind*0.5;
end

s.Processed.PhotoRaster = photoRaster;
s.Processed.VelRaster = velRaster;
s.Processed.PhotoRasterZ = photoRasterZ;
s.Processed.BaselineMean = baselineMean;
s.Processed.BaselineSTD = baselineSTD;


%now I need to store different means

photoAverages = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,numFreqs,numDBs);
velAverages = zeros(rasterVelWindow(2) - rasterVelWindow(1) + 1,numFreqs,numDBs);
photoRasterStore = cell(numFreqs,numDBs);
velRasterStore = cell(numFreqs,numDBs);
for ind = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(ind),3));
    for j = 1:numDBs
        %find all trials of the particular setting
        targetFinder = find(trialMatrix(:,2) == uniqueFreqs(ind) & trialMatrix(:,3) == subUniqueDB(j));
        %pull and average these traces
        tempHolder = photoRasterZ(:,targetFinder);
        photoRasterStore{ind,j} = tempHolder;
        tempHolder = mean(tempHolder');
        photoAverages(:,ind,j) = tempHolder;
        %pull velocity traces and do the same
        tempVel = velRaster(:,targetFinder);
        velRasterStore{ind,j} = tempVel;
        tempVel = mean(tempVel');
        velAverages(:,ind,j) = tempVel;
        
    end
end

s.Processed.PhotoAverages = photoAverages;
s.Processed.PhotoStore = photoRasterStore;
s.Processed.VelAverages = velAverages;
s.Processed.VelStore = velRasterStore;

%lets make rasters with the peak detection rise times (based on Chris'
%recommendation)

%make basic rasters
[riseRasters] = functionBasicRaster((targetPeaks(:,5)),traceMBED,rasterWindow);
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

%now lets find the indices for various things to sort out by amplitude

%first, lets create an alternate trialMatrix with standardized DBs
standardDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(1),3));
altMatrix = trialMatrix;
for ind = 2:numFreqs
    findDBs = unique(altMatrix(altMatrix(:,2) == uniqueFreqs(ind),3));
    for j = 1:length(findDBs)
        altMatrix(altMatrix(:,2) == uniqueFreqs(ind) & altMatrix(:,3) == findDBs(j),3) = standardDB(j);
    end
end

%now using the altMatrix, find DB related stuff!
dbSort = zeros(length(altMatrix),numDBs);
for ind = 1:numDBs
    sortingCounter = 1;
    sortIndex = zeros(length(altMatrix),1);
    for j = 1:numFreqs
        
        sortingFinder = find(altMatrix(:,2) == uniqueFreqs(j) & altMatrix(:,3) == uniqueDBs(ind));
        sortIndex(sortingFinder) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
    dbSort(:,ind) = sortIndex;
    sortIndex = [];
end

for i = 1:numDBs
    riseRasters(:,3+i) = dbSort(riseRasters(:,2),i);
end

for i = 1:numFreqs
    for j = 1:numDBs
        tempData = s.Processed.PhotoStore{i,j};
        baseVal = min(tempData(130:145,:));
        peakVal = max(tempData(145:192,:));
        baseMean = mean(tempData(130:145,:));
        peakMean = mean(tempData(145:192,:));
        peakMag = peakVal - baseVal;
        peakInt = peakMean - baseMean;
        magStore(:,i,j) = peakMag;
        intStore(:,i,j) = peakInt;
    end
end


s.Processed.RiseRaster = riseRasters;
s.Processed.MagStore = magStore;
s.Processed.IntStore = intStore;
s.SoundData.AltMatrix = altMatrix;
s.SoundData.DBSort = dbSort;



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
if isfield(soundData,'WhiteNoise')
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
rasterAxis(:,2) = [0:size(photoRasterZ,1)/8:size(photoRasterZ,1)];
rasterAxis(1,2) = 1;


%create a set of points for the velRaster axis
velRasterAxis = zeros(rasterVelWindow(2)/10-rasterVelWindow(1)/10+1,2);
velRasterAxis(:,1) = [rasterVelWindow(1)/10:rasterVelWindow(2)/10];
%find a second in photometry sample time
velRasterAxis(:,2) = [0:length(velAverages)/8:length(velAverages)];
velRasterAxis(1,2) = 1;


%also create labels for the rasters. first column is frequencies, second is
%where lines are drawn, third is for placement of labels
rasterLabels = zeros(numFreqs,3);
rasterLabels(:,1) = uniqueFreqs;
rasterLabels(:,2) = [numDBs*soundData.ToneRepetitions:numDBs*soundData.ToneRepetitions:length(soundData.TrialMatrix)];
rasterLabels(:,3) = [numDBs*soundData.ToneRepetitions/2:numDBs*soundData.ToneRepetitions:length(soundData.TrialMatrix)-numDBs*soundData.ToneRepetitions/2];

dbRasterLabels = zeros(numFreqs,3);
dbRasterLabels(:,1) = uniqueFreqs;
dbRasterLabels(:,2) = [soundData.ToneRepetitions:soundData.ToneRepetitions:length(soundData.TrialMatrix)/numDBs];
dbRasterLabels(:,3) = [soundData.ToneRepetitions/2:soundData.ToneRepetitions:length(soundData.TrialMatrix)/numDBs-soundData.ToneRepetitions/2];


rasterTicks = cell(numFreqs,1);
for ind = 1:numFreqs
    rasterTicks{ind} = num2str(rasterLabels(ind,1));
end


%MAKE THE FIGURE
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])


imagescLim = [min(min(min(photoAverages))),max(max(max(photoAverages)))];

velMax = max(max(max(velAverages)));
velMin = min(min(min(velAverages)));

for i = 1:numDBs
    %plot heatmap of all responses
    subplot(numDBs,4,1+numDBs*(i-1))
    imagesc(squeeze(photoAverages(:,:,i)'),imagescLim)
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    if i == 1
        title({fileName;strcat(num2str(uniqueDBs(i)),'DB')})
    else
        title({strcat(num2str(uniqueDBs(i)),'DB')})
    end
    
    %plot peak detection rasters
    subplot(numDBs,4,2+numDBs*(i-1))
    hold on
    plot(riseRasters(:,1),riseRasters(:,3+i),'k.')

    for ind = 1:numFreqs
        plot([rasterWindow(1) rasterWindow(2)],[dbRasterLabels(ind,2) dbRasterLabels(ind,2)],'g')
    end
    plot([0 0],[1 length(trialMatrix)],'b')
    plot([soundData.ToneDuration soundData.ToneDuration],[1 length(trialMatrix)],'b')
    title(strcat(num2str(uniqueDBs(i)),'RastersByFreq'))
    set(gca,'YTick',dbRasterLabels(:,3))
    set(gca,'YTickLabel',rasterTicks)
    ylim([2 s.SoundData.ToneRepetitions*numFreqs])
    xlim(rasterWindow)
    set(gca,'Ydir','reverse')
    
    %plot out velocity changes
    subplot(numDBs,4,3+numDBs*(i-1))
    imagesc(squeeze(velAverages(:,:,i))',[velMin velMax])
    colormap('parula')
    set(gca,'XTick',velRasterAxis(:,2));
    set(gca,'XTickLabel',velRasterAxis(:,1));
    title(strcat('Velocity Heatmap',num2str(uniqueDBs(i)),'DB'))
    
    %plot out individual magnitudes
    subplot(numDBs,4,4+numDBs*(i-1))
    hold on
    plot(squeeze(magStore(:,:,i)))
    meanFind = mean(squeeze(magStore(:,:,i))');
    plot(meanFind,'k','LineWidth',2)
    title('Individual Response Magnitudes')
    ylim([0 max(max(squeeze(magStore(:,:,i))))])
    xlim([1,numReps])
    
end
% 
% %plot overall photometry trace and locomotion trace
% subplot(3,4,1)
% imagesc(squeeze(photoAverages(:,:,1)'),imagescLim)
% colormap('parula')
% colorbar
% set(gca,'XTick',rasterAxis(:,2));
% set(gca,'XTickLabel',rasterAxis(:,1));
% set(gca,'YTick',octaveRange(:,2));
% set(gca,'YTickLabel',octaveRange(:,1));
% title({fileName;'60DB'})
% 
% subplot(3,4,5)
% imagesc(squeeze(photoAverages(:,:,2)'),imagescLim)
% colormap('parula')
% colorbar
% set(gca,'XTick',rasterAxis(:,2));
% set(gca,'XTickLabel',rasterAxis(:,1));
% set(gca,'YTick',octaveRange(:,2));
% set(gca,'YTickLabel',octaveRange(:,1));
% title('80DB')
% 
% subplot(3,4,9)
% imagesc(squeeze(photoAverages(:,:,3)'),imagescLim)
% colormap('parula')
% colorbar
% set(gca,'XTick',rasterAxis(:,2));
% set(gca,'XTickLabel',rasterAxis(:,1));
% set(gca,'YTick',octaveRange(:,2));
% set(gca,'YTickLabel',octaveRange(:,1));
% title('100DB')
% 
% %plot out rasters
% subplot(3,4,2)
% hold on
% plot(riseRasters(:,1),riseRasters(:,4),'k.')
% 
% for ind = 1:numFreqs
%     plot([rasterWindow(1) rasterWindow(2)],[dbRasterLabels(ind,2) dbRasterLabels(ind,2)],'g')
% end
% plot([0 0],[1 length(trialMatrix)],'b')
% plot([soundData.ToneDuration soundData.ToneDuration],[1 length(trialMatrix)],'b')
% title('60DB Rasters Organized by Freq')
% set(gca,'YTick',dbRasterLabels(:,3))
% set(gca,'YTickLabel',rasterTicks)
% ylim([2 s.SoundData.ToneRepetitions*numFreqs])
% xlim(rasterWindow)
% set(gca,'Ydir','reverse')
% 
% subplot(3,4,6)
% hold on
% plot(riseRasters(:,1),riseRasters(:,5),'k.')
% 
% for ind = 1:numFreqs
%     plot([rasterWindow(1) rasterWindow(2)],[dbRasterLabels(ind,2) dbRasterLabels(ind,2)],'g')
% end
% plot([0 0],[1 length(trialMatrix)],'b')
% plot([soundData.ToneDuration soundData.ToneDuration],[1 length(trialMatrix)],'b')
% title('80DB Rasters Organized by Freq')
% set(gca,'YTick',dbRasterLabels(:,3))
% set(gca,'YTickLabel',rasterTicks)
% ylim([2 s.SoundData.ToneRepetitions*numFreqs])
% xlim(rasterWindow)
% set(gca,'Ydir','reverse')
% 
% subplot(3,4,10)
% hold on
% plot(riseRasters(:,1),riseRasters(:,6),'k.')
% 
% for ind = 1:numFreqs
%     plot([rasterWindow(1) rasterWindow(2)],[dbRasterLabels(ind,2) dbRasterLabels(ind,2)],'g')
% end
% plot([0 0],[1 length(trialMatrix)],'b')
% plot([soundData.ToneDuration soundData.ToneDuration],[1 length(trialMatrix)],'b')
% title('100DB Rasters Organized by Freq')
% set(gca,'YTick',dbRasterLabels(:,3))
% set(gca,'YTickLabel',rasterTicks)
% ylim([2 s.SoundData.ToneRepetitions*numFreqs])
% xlim(rasterWindow)
% set(gca,'Ydir','reverse')
% %plot out velocities
% 
% %find min and max overall
% 
% 
% subplot(3,4,3)
% imagesc(squeeze(velAverages(:,:,1))',[velMin velMax])
% colormap('parula')
% set(gca,'XTick',velRasterAxis(:,2));
% set(gca,'XTickLabel',velRasterAxis(:,1));
% title('Velocity Heatmap 60DB')
% 
% subplot(3,4,7)
% imagesc(squeeze(velAverages(:,:,2))',[velMin velMax])
% colormap('parula')
% set(gca,'XTick',velRasterAxis(:,2));
% set(gca,'XTickLabel',velRasterAxis(:,1));
% title('Velocity Heatmap 80DB')
% 
% subplot(3,4,11)
% imagesc(squeeze(velAverages(:,:,3))',[velMin velMax])
% colormap('parula')
% set(gca,'XTick',velRasterAxis(:,2));
% set(gca,'XTickLabel',velRasterAxis(:,1));
% title('Velocity Heatmap 100DB')

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


diary off
end


