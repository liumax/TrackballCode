%This is code that will read and link MBED inputs with photometry, while
%recording from my big rig. 

% function [s] = analysisPhotometryTuning(fileName);
fileName = '180105_ML171117A_Tuning';
% fileName = '180105_ML171117B_Tuning';
% fileName = '180105_ML171119A_Tuning';
% fileName = '180105_ML171119B_Tuning';

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

%% Extracts Sound Data from soundFile, including freq, db, reps.
soundName = strcat(fileName,'Sound.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;

%check for presence of white noise
try
    whiteStatus = soundData.WhiteNoise;
catch
    whiteStatus = 0;
end

%pull important sound and trial information
uniqueFreqs = unique(soundData.Frequencies);
uniqueDBs = unique(soundData.TrialMatrix(soundData.TrialMatrix(:,2)==uniqueFreqs(1),3));
uniqueDBSteps = unique(soundData.TrialMatrix(soundData.TrialMatrix(:,2)==uniqueFreqs(1),3));
numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);
trialMatrix = soundData.TrialMatrix;
expectedITI = soundData.Delays;

%store these in structured array
s.SoundData.UniqueFrequencies = uniqueFreqs;
s.SoundData.UniqueDBs = uniqueDBs;
s.SoundData.NumFreqs = numFreqs;
s.SoundData.NumDBs = numDBs;

toneDur = soundData.ToneDuration;
toneReps = soundData.ToneRepetitions;
totalTrialNum = length(soundData.Frequencies);


% Generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for i =1:numFreqs
    freqNameHolder{i} = num2str(uniqueFreqs(i));
end

if whiteStatus ==0
    % Finds the number of octaves, makes array of octave steps. This will be used for imagesc graphing applications
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
elseif whiteStatus == 1;
    % Finds the number of octaves, makes array of octave steps. This will be used for imagesc graphing applications
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
    octaveRange(2:end+1,:) = octaveRange(1:end,:);
    octaveRange(1,1) = 0;
    octaveRange(1,2) = 1;
end

% Does the same for dBs. 
if length(uniqueDBs) == 1
    dbSteps = 1;
    totalDBs = 1;
    dbRange = [100,1];
else
    dbSteps = uniqueDBSteps(2) - uniqueDBSteps(1);
    totalDBs = (uniqueDBSteps(end) - uniqueDBSteps(1))/dbSteps;
    dbRange = zeros(totalDBs + 1,2);
    dbRange(:,1) = uniqueDBSteps(1):dbSteps:uniqueDBSteps(end);
    for i = 1:size(dbRange,1)
        dbRange(i,2) = find(uniqueDBSteps == dbRange(i,1));
    end
end


%stores into s.
s.Parameters.OctaveRange = octaveRange;
s.Parameters.DBRange = dbRange;

disp('Sound Data Extracted')
%% Next, lets pull the MBED stuff
disp('Loading MBED...')
[trialStates, portStates, trialParams] = maxTrialVariablesLickingTask(fileName);
disp('MBED Data Loaded!')
%now extract ports 1 (tone) and port 7 (inscopix sync)

mbedTone = [portStates.tStamps',portStates.inStates(:,1)];
mbedSync = [portStates.tStamps',portStates.inStates(:,7)];

[mbedTone] = functionSignalDuplicateElim(mbedTone,2);
[mbedSync] = functionSignalDuplicateElim(mbedSync,2);

s.MBED.ToneDelivery = mbedTone;
s.MBED.Sync = mbedSync;

frameTimes = mbedSync(mbedSync(:,2) == 1,1);
frameRate = 1/(mean(diff(frameTimes))/1000);

%repair frameTimes!
meanITI = mean(diff(frameTimes));
whileTrig = 0;
counter = 1;
frameDiffs = diff(frameTimes);
while whileTrig == 0
    if counter == length(frameDiffs)
        break
    end
    if frameDiffs(counter) > meanITI*1.3
        disp(strcat('Missed Frame Found:',num2str(counter)))
        frameTimes(counter+2:end+1) = frameTimes(counter+1:end);
        frameTimes(counter+1) = frameTimes(counter) + meanITI;
        frameDiffs = diff(frameTimes);
        disp('Missing Frame Repaired')
    end
    counter = counter + 1;
end

toneTimes = mbedTone(mbedTone(:,2) == 1,1);
%CHECK FOR ISSUES BETWEEN ACTUAL SOUND DATA AND RECEIVED TONES
if length(toneTimes) == length(soundData.Delays)
    disp('Tone Times and Tuning Delay Number Matched')
else
    disp('Mismatch in Number of Tone Times and Tuning Delays')
    lengthDiff = length(toneTimes) - length(soundData.Delays);
    if lengthDiff>0
        %first, check to see if there are big ITIs
        diffTone = diff(toneTimes);
        maxDiff = max(expectedITI)*1000*1.1;
        diffFind = find(diffTone>maxDiff);
        toneTimes(diffFind) = [];
        disp('Removed Extra TTL')
        disp(strcat('MBED TTLs:',num2str(length(toneTimes))))
        disp(strcat('File TTLs:',num2str(length(soundData.Delays))))
        figure
        hold on
        plot(diff(toneTimes)/1000)
        plot(expectedITI,'r')
        title('Tone ITIs MBED Blue Expected R')
        xlabel('Tone Number')
        ylabel('Seconds')
    else
        error('Unplanned Error')
    end
end


%Now pull locomotor data, only if no TMP file
[locoData] = functionLocoTmp(fileName,portStates,locoTimeStep);
s.Locomotion = locoData;
disp('Locomotion Data Loaded')

%% Now pull up CSV file!
trace = csvread(strcat(fileName,'.csv'),1,1);

%% now lets make basic raster!

if length(trace) == length(frameTimes)
    disp('Photometry Trace and MBED Frames Equal')
elseif length(trace) ~= length(frameTimes)
    error('Mismatch in Photometry Trace and MBED Frames')
end

%now we need to convert things to inscopix frames
alignTimes = zeros(length(toneTimes),1);
for i = 1:length(alignTimes)
    finder = find(frameTimes - toneTimes(i) >=0,1,'first');
    alignTimes(i) = finder;
end

rasterPhotWindow = round(rasterWindow*frameRate);

for i = 1:length(alignTimes)
    photoRaster(:,i) = trace(alignTimes(i) + rasterPhotWindow(1):alignTimes(i) + rasterPhotWindow(2));
end

[photoRasterZ,baselineMean,baselineSTD] = functionZScore(photoRaster,rasterWindow(1)*-1,length(alignTimes));

rasterVelWindow = round(rasterWindow/locoTimeStep);
velRaster = zeros(rasterVelWindow(2) - rasterVelWindow(1) + 1,length(alignTimes));
for ind = 1:length(toneTimes)
    %find the time in the velocity trace
    alignTime = toneTimes(ind)/1000;
    velPoint = find(locoData.Velocity(:,1) - alignTime > 0,1,'first');
    if velPoint + rasterVelWindow(2) < length(locoData.Velocity(:,1)) & velPoint + rasterVelWindow(1) >= 1
        velRaster(:,ind) = locoData.Velocity(velPoint + rasterVelWindow(1):velPoint + rasterVelWindow(2),2); 
    else
        disp(ind)
        disp('Vel Rasters: Not matched up')
    end
end

%make simple plot
simplePlot = photoRaster;
for ind = 1:length(toneTimes)
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

% for i = 1:numFreqs
%     for j = 1:numDBs
%         tempData = s.Processed.PhotoStore{i,j};
%         baseVal = min(tempData(130:145,:));
%         peakVal = max(tempData(145:192,:));
%         baseMean = mean(tempData(130:145,:));
%         peakMean = mean(tempData(145:192,:));
%         peakMag = peakVal - baseVal;
%         peakInt = peakMean - baseMean;
%         magStore(:,i,j) = peakMag;
%         intStore(:,i,j) = peakInt;
%     end
% end

% 
% 
% s.Processed.MagStore = magStore;
% s.Processed.IntStore = intStore;
s.SoundData.AltMatrix = altMatrix;
s.SoundData.DBSort = dbSort;


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
%     
%     %plot peak detection rasters
%     subplot(numDBs,4,2+numDBs*(i-1))
%     hold on
%     plot(riseRasters(:,1),riseRasters(:,3+i),'k.')
% 
%     for ind = 1:numFreqs
%         plot([rasterWindow(1) rasterWindow(2)],[dbRasterLabels(ind,2) dbRasterLabels(ind,2)],'g')
%     end
%     plot([0 0],[1 length(trialMatrix)],'b')
%     plot([soundData.ToneDuration soundData.ToneDuration],[1 length(trialMatrix)],'b')
%     title(strcat(num2str(uniqueDBs(i)),'RastersByFreq'))
%     set(gca,'YTick',dbRasterLabels(:,3))
%     set(gca,'YTickLabel',rasterTicks)
%     ylim([2 s.SoundData.ToneRepetitions*numFreqs])
%     xlim(rasterWindow)
%     set(gca,'Ydir','reverse')
    
    %plot out velocity changes
    subplot(numDBs,4,3+numDBs*(i-1))
    imagesc(squeeze(velAverages(:,:,i))',[velMin velMax])
    colormap('parula')
    set(gca,'XTick',velRasterAxis(:,2));
    set(gca,'XTickLabel',velRasterAxis(:,1));
    title(strcat('Velocity Heatmap',num2str(uniqueDBs(i)),'DB'))
    
%     %plot out individual magnitudes
%     subplot(numDBs,4,4+numDBs*(i-1))
%     hold on
%     plot(squeeze(magStore(:,:,i)))
%     meanFind = mean(squeeze(magStore(:,:,i))');
%     plot(meanFind,'k','LineWidth',2)
%     title('Individual Response Magnitudes')
%     ylim([0 max(max(squeeze(magStore(:,:,i))))])
%     xlim([1,numReps])
    
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
% end


