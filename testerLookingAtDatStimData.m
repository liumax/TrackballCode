











%this is meant to be master code to run analysis code!
%Establishes the home folder as the base
masterFolder = pwd;
masterDir = dir;
masterDir = {masterDir.name};
%this code here is to remove folders reported in dir like "." or ".."
masterIndex = strfind(masterDir,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
%masterFolders is a listing of all folders that I want to target for batch
%analysis
masterFolders = masterDir(masterIndex);
numFolders = size(masterFolders,2);
for masterCount = 1:numFolders
%generates text that signifies path to the target folder
targetPath = strcat(masterFolder,'\',masterFolders{masterCount});
cd(targetPath)
%what lists matlab files in the folder. can extract based on file type.
test = what;
test = test.mat;
test = test{1};
load(test);
nameFind = strfind(test,'DatStimAnalysis');
newName =strcat('x',test(1:nameFind-1));
fullData.(newName) = s;
cd(masterFolder)
end


%pull names of recordings out
fieldNames = fields(fullData);
%pull individual components of names out
fieldData = cell(length(fieldNames),5);
for i = 1:length(fieldNames)
    %pull name as string
    nameHold = fieldNames{i};
    spaceFind = strfind(nameHold,'_');
    %store date
    fieldData{i,1} = nameHold(2:spaceFind(1)-1);
    %store animal name
    fieldData{i,2} = nameHold(spaceFind(1)+1:spaceFind(2)-1);
    %store hemisphere
    fieldData{i,3} = nameHold(spaceFind(2)+1:spaceFind(3)-1);
    %store depth
    fieldData{i,4} = nameHold(spaceFind(3)+1:spaceFind(4)-1);
    %store type
    fieldData{i,5} = nameHold(spaceFind(4)+1:end);
end




store.FileName = cell(10,1);
store.UnitName = cell(10,1);
store.BaselineRate = zeros(10,1);
store.BaselineSTD = zeros(10,1);
store.PostZeroBin = zeros(10,1);
store.PostZeroSTD = zeros(10,1);
store.PostBigBin = zeros(10,1);
store.TrueAUC = zeros(10,1);
store.ShuffleAUC = zeros(10,2);
store.Waves = zeros(40,4,10);
placeHolder = 1;
for i = 1:length(fieldNames)
    %load data
    testData = fullData.(fieldNames{i});
    desigNames = testData.DesignationName;
    %pull tuning data
    %go into individual units to pull data
    for j = 1:length(desigNames)
        store.FileName{placeHolder - 1 + j} = fieldNames{i};
        store.UnitName{placeHolder - 1 + j} = desigNames{j};
        store.StartPref(placeHolder - 1 + j) = testData.RotaryData.StartPreference;
        store.EndPref(placeHolder - 1 + j) = testData.RotaryData.EndPreference;
        store.PrefRange((placeHolder - 1 + j),1) = prctile(testData.RotaryData.PreferenceDistribution,0.5);
        store.PrefRange((placeHolder - 1 + j),2) = prctile(testData.RotaryData.PreferenceDistribution,99.5);
        store.BaselineRate(placeHolder - 1 + j) = testData.(desigNames{j}).AverageRate;
        store.BaselineSTD(placeHolder - 1 + j) = testData.(desigNames{j}).AverageSTD;
        store.PostZeroBin(placeHolder - 1 + j) = testData.(desigNames{j}).AllHistograms(81);
        store.PostZeroSTD(placeHolder - 1 + j) = testData.(desigNames{j}).HistogramStandardDeviation(81);
        store.PostBigBin(placeHolder - 1 + j) = mean(testData.(desigNames{j}).AllHistograms(81:101));
        store.FullHist(:,placeHolder - 1 + j) = testData.(desigNames{j}).AllHistograms;
        if isfield(testData.(desigNames{j}),'TrueAUC')
        store.TrueAUC(placeHolder - 1 + j) = testData.(desigNames{j}).TrueAUC;
        store.ShuffleAUC(placeHolder - 1 + j,1) = prctile(testData.(desigNames{j}).ShuffleAUC,0.5);
        store.ShuffleAUC(placeHolder - 1 + j,2) = prctile(testData.(desigNames{j}).ShuffleAUC,99.5);
        end
        %store waveform
        store.Waves(:,:,placeHolder - 1 + j) = testData.(desigNames{j}).AverageWaveForms;
    end
    placeHolder = placeHolder + length(desigNames);
end


%lets plot some basics

hFig = figure
hist(store.BaselineRate,100);
title('Distribution of Baseline Firing Rates For Terminal Stimulation Experiments')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'BaselineHistogramTerminalStim','-dpdf','-r0')


for i = 1:placeHolder - 1
    if store.TrueAUC(i) > store.ShuffleAUC(i,2) | store.TrueAUC(i) < store.ShuffleAUC(i,1)
        findSigAUC(i) = 1;
    else
        findSigAUC(i) = 0;
    end
end
findSigAUC = find(findSigAUC == 1);

hFig = figure
plot(store.TrueAUC,store.BaselineRate,'k.')
hold on
plot(store.TrueAUC(findSigAUC),store.BaselineRate(findSigAUC),'ro')
title({'AUC Values vs Baseline Rate';'Red Circles are Significant by 99% CI';strcat(num2str(length(findSigAUC)/(placeHolder-1)),'% Units Significant')})
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'BaselineVsAUCTerminal','-dpdf','-r0')

%look at post big bin and post first bin
for i = 1:placeHolder-1
    postZeroZ(i) = (store.PostZeroBin(i)-store.BaselineRate(i))/store.BaselineSTD(i);
    postBigZ(i) = (store.PostBigBin(i)-store.BaselineRate(i))/store.BaselineSTD(i);
end

hFig = figure
hist(postZeroZ,100)
title('Histogram of Z-Score of First 50ms bin following onset of laser stimulation')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'postZeroZTerminal','-dpdf','-r0')

hFig = figure
hist(postBigZ,100)
title('Histogram of Z-Score of First 1s bin following onset of laser stimulation')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'postBigZTerminal','-dpdf','-r0')

%combine all units
hFig = figure
plot(mean(store.FullHist'),'LineWidth',2)
hold on
plot(mean(store.FullHist')+std(store.FullHist')/sqrt(placeHolder-1))
plot(mean(store.FullHist')-std(store.FullHist')/sqrt(placeHolder-1))
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Averaged Firing from All Units, Aligned to Laser')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'AverageAllUnitsTerminal','-dpdf','-r0')


%separate by recording
for i = 1:length(fieldNames)
    test = strfind(store.FileName,fieldNames{i});
    fieldStore = find(not(cellfun('isempty', test)));
    fieldTicker(fieldStore) = i;
end

%how about histograms of responses?
for i = 1:placeHolder - 1
    normHist(:,i) = store.FullHist(:,i);
    zHist(:,i) = (normHist(:,i)-store.BaselineRate(i))/store.BaselineSTD(i);
    smoothZ(:,i) = smooth(zHist(:,i),9);
    histMin = min(normHist(:,i));
    histMax = max(normHist(:,i));
    normHist(:,i) = (normHist(:,i)-histMin)/(histMax-histMin); 
    smoothNorm(:,i) = smooth(normHist(:,i),9);
end

%for plotting separated by field/experiment
% for i = 1:length(fieldNames)
%     testData = normHist(:,fieldTicker == i);
%     testBaseline = store.BaselineRate(fieldTicker==i);
%     [B,I] = sort(testBaseline);
%     figure
%     imagesc(testData(:,I)')
% end

%order from highest baseline 
[B,I] = sort(store.BaselineRate);
newZHist = zHist(:,I);
newNormHist = normHist(:,I);

hFig = figure
imagesc(newNormHist')
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Normalized Responses Sorted by Baseline Firing Rate')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapNormHistTerminal','-djpeg','-r0')

hFig = figure
imagesc(newZHist',[-2,2])
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Z-Score Responses Sorted by Baseline Firing Rate')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapZHistTerminal','-djpeg','-r0')


hFig = figure
imagesc(smoothNorm(:,I)')
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Normalized Responses Sorted by Baseline Firing Rate')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapSmoothNormHistTerminal','-djpeg','-r0')

hFig = figure
imagesc(smoothZ(:,I)',[-2,2])
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Z-Score Responses Sorted by Baseline Firing Rate')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapSmoothZHistTerminal','-djpeg','-r0')
% in zscored see some pretty long lasting inhibition in the fast spiking
% cells. Could just be an artifact of the laser stimulation?

%okay, now lets get at the running. DO we see preference for running or
%stopping?

for i = 1:placeHolder - 1
    if store.StartPref(i) > store.PrefRange(i,2)
        findSigStart(i) = 1;
    else
        findSigStart(i) = 0;
    end
end
findSigStart = find(findSigStart == 1);


for i = 1:placeHolder - 1
    if store.EndPref(i) > store.PrefRange(i,2)
        findSigEnd(i) = 1;
    else
        findSigEnd(i) = 0;
    end
end
findSigEnd = find(findSigEnd == 1);

%looks like no changes for locomotion starts or ends. What about looking at
%average locomotion trace?
for i = 1:length(fieldNames)
    testData = fullData.(fieldNames{i});
    %lets make average locomotion trace!
    rasterWindow = testData.Parameters.RasterWindow;
    %pull laser times
    laserTimes = testData.LaserData.LaserStartTimes;
    laserTimeDiff = min(diff(laserTimes));
    laserTimeStarts = [1;find(diff(laserTimes)>laserTimeDiff*2)+1];
    
    for j = 1:length(laserTimeStarts)
        targetTime = laserTimes(laserTimeStarts(j));
        %find start and end
        velStart = find(testData.RotaryData.Velocity(:,1)<targetTime+rasterWindow(1),1,'last');
        velEnd = velStart+((rasterWindow(2)-rasterWindow(1))/0.01);
        velSnip = testData.RotaryData.Velocity(velStart:velEnd,2);
        velHold(:,j) = velSnip;
    end
    velRaster{i} = velHold;
    velHold = [];
    laserReps(i) = length(laserTimeStarts);
    velMean(:,i) = mean(velRaster{i}');
    velSTD(:,i) = std(velRaster{i}');
    velSTE(:,i) = velSTD(:,i)/sqrt(length(laserTimeStarts));
end

%dont see much in the average velocity traces either lol


hFig = figure
hold on
for i = 1:size(velMean,2)
    plot(smooth(velMean(:,i)-mean(velMean(:,i))+(2*(i-1)),50))
end
title('Average Velocity Traces')
set(gca,'XTick',[0:100:1000])
set(gca,'XTickLabel',[-4:1:6])
xlim([0,1000])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'AverageVelocityTerminal','-dpdf','-r0')

%lets try normalizing
for i = 1:size(velMean,2)
    normVel(:,i) = velMean(:,i);
    zVel(:,i) = (normVel(:,i) - mean(normVel(:,i)))/std(normVel(:,i));
    zVel(:,i) = smooth(zVel(:,i),50);
    velMax = max(normVel(:,i));
    velMin = min(normVel(:,i));
    normVel(:,i) = (normVel(:,i)-velMin)/(velMax-velMin);
    normVel(:,i) = smooth(normVel(:,i),50);
end

hFig = figure
imagesc(flipud(zVel'),[-2 2])
colorbar
title('Z-Scored Velocity')
set(gca,'XTick',[0:100:1000])
set(gca,'XTickLabel',[-4:1:6])
xlim([0,1000])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'ZVelocitySmoothTerminal','-dpdf','-r0')











